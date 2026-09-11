function results = lcpm_repeated_cv(x, y, covariates, cfg, varargin)
%LCPM_REPEATED_CV Repeated k-fold cross-validation of the LCPM classifier.
% Inputs:
%   x          - n_edges x n_subjects matrix of Fisher-z connectivity
%                values (from lcpm_prepare_edges)
%   y          - n_subjects x 1 binary labels (0/1)
%   covariates - n_subjects x n_covariates matrix
%   cfg        - parameter struct from config.m
%
% Optional name-value pairs:
%   'RandInds' - n_subjects x cfg.n_iterations matrix of predefined random
%                permutations. When supplied, the exact same data splits
%                are reused in every iteration (used by the permutation
%                test and the virtual lesion analysis).
%   'EdgeMask' - logical vector selecting a subset of edges. When
%                supplied, the classifier is trained on the selected edges
%                only and edge-selection frequencies are not collected
%                (used by the virtual lesion analysis).
%
% Pipeline (per fold, training data only):
%   1. Feature selection: point-biserial (Pearson) correlation between
%      each edge and the label; edges with p < cfg.p_threshold are kept and
%      split into positively / negatively correlated sets.
%   2. Composite feature: sum of positive edges minus sum of negative
%      edges.
%   3. Weighted logistic regression on the composite feature plus
%      covariates, standardized on the training folds, with weights
%      inversely proportional to class size.
%   4. Evaluation on the held-out fold: AUC (perfcurve) and
%      accuracy/sensitivity/specificity at the fixed probability cutoff
%      cfg.class_threshold.
%
% Results:
%   results - struct with per-iteration metrics (auc/acc/sen/spec/edges),
%             edge-selection frequencies, the random splits used, and the
%             mean/std summary over valid iterations.

p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addRequired(p, 'covariates', @isnumeric);
addRequired(p, 'cfg', @isstruct);
addParameter(p, 'RandInds', [], @isnumeric);
addParameter(p, 'EdgeMask', [], @(v) isempty(v) || isvector(v));
parse(p, x, y, covariates, cfg, varargin{:});

x          = p.Results.x;
y          = double(p.Results.y(:));
covariates = p.Results.covariates;
cfg        = p.Results.cfg;
rand_inds_all = p.Results.RandInds;
edge_mask     = p.Results.EdgeMask;

if ~ismatrix(x)
    error('x must be a 2-D edge x subject matrix (use lcpm_prepare_edges for 3-D connectivity arrays).');
end
[n_edges_full, n_subjects] = size(x);

if ~isempty(edge_mask)
    edge_mask = logical(edge_mask(:));
    if numel(edge_mask) ~= n_edges_full
        error('EdgeMask has %d entries but x has %d edges.', numel(edge_mask), n_edges_full);
    end
    x = x(edge_mask, :);
end

n_iterations = cfg.n_iterations;
n_folds      = cfg.n_folds;

% Determine the random splits. By default iteration i uses rng(i) so the
% splits are identical to the original analysis scripts; supplying
% 'RandInds' reuses predefined splits instead.
if isempty(rand_inds_all)
    rand_inds_all = zeros(n_subjects, n_iterations);
    for i = 1:n_iterations
        rng(i, 'twister');
        rand_inds_all(:, i) = randperm(n_subjects);
    end
else
    if ~isequal(size(rand_inds_all), [n_subjects, n_iterations])
        error('RandInds must be of size n_subjects x cfg.n_iterations (%d x %d).', ...
            n_subjects, n_iterations);
    end
end

% Edge-selection frequencies are only meaningful for the full edge set.
collect_counts = isempty(edge_mask);
if collect_counts
    edge_selected_count     = zeros(n_edges_full, 1);
    pos_edge_selected_count = zeros(n_edges_full, 1);
    neg_edge_selected_count = zeros(n_edges_full, 1);
else
    edge_selected_count     = [];
    pos_edge_selected_count = [];
    neg_edge_selected_count = [];
end

auc_per_iter  = nan(n_iterations, 1);
acc_per_iter  = nan(n_iterations, 1);
sen_per_iter  = nan(n_iterations, 1);
spec_per_iter = nan(n_iterations, 1);
edges_per_iter = nan(n_iterations, 1);

base_size = floor(n_subjects / n_folds);
remainder = n_subjects - base_size * n_folds;

for iter = 1:n_iterations
    randinds  = rand_inds_all(:, iter);
    start_idx = 1;

    auc_folds   = nan(n_folds, 1);
    acc_folds   = nan(n_folds, 1);
    sen_folds   = nan(n_folds, 1);
    spec_folds  = nan(n_folds, 1);
    edges_folds = nan(n_folds, 1);

    for fold = 1:n_folds
        test_size = base_size + (fold <= remainder);
        test_pos  = start_idx : (start_idx + test_size - 1);
        start_idx = start_idx + test_size;
        train_pos = [1:test_pos(1) - 1, test_pos(end) + 1:n_subjects];

        test_inds  = randinds(test_pos);
        train_inds = randinds(train_pos);

        X_train = x(:, train_inds).';
        y_train = y(train_inds);
        cov_train = covariates(train_inds, :);
        X_test  = x(:, test_inds).';
        y_test  = y(test_inds);
        cov_test = covariates(test_inds, :);

        % Feature selection on the training folds only.
        [rho, pval] = corr(X_train, y_train, 'Type', 'Pearson');
        pos_idx = find(pval < cfg.p_threshold & rho > 0);
        neg_idx = find(pval < cfg.p_threshold & rho < 0);

        edges_folds(fold) = numel(pos_idx) + numel(neg_idx);

        if collect_counts
            edge_selected_count(pos_idx)     = edge_selected_count(pos_idx) + 1;
            edge_selected_count(neg_idx)     = edge_selected_count(neg_idx) + 1;
            pos_edge_selected_count(pos_idx) = pos_edge_selected_count(pos_idx) + 1;
            neg_edge_selected_count(neg_idx) = neg_edge_selected_count(neg_idx) + 1;
        end

        if isempty(pos_idx) && isempty(neg_idx)
            continue;
        end

        % Composite feature: sum of positive edges minus sum of negative edges.
        if isempty(pos_idx)
            pos_sum_train = zeros(size(X_train, 1), 1);
            pos_sum_test  = zeros(size(X_test, 1), 1);
        else
            pos_sum_train = sum(X_train(:, pos_idx), 2);
            pos_sum_test  = sum(X_test(:, pos_idx), 2);
        end
        if isempty(neg_idx)
            neg_sum_train = zeros(size(X_train, 1), 1);
            neg_sum_test  = zeros(size(X_test, 1), 1);
        else
            neg_sum_train = sum(X_train(:, neg_idx), 2);
            neg_sum_test  = sum(X_test(:, neg_idx), 2);
        end
        feat_train = pos_sum_train - neg_sum_train;
        feat_test  = pos_sum_test  - neg_sum_test;

        % Standardize the composite feature and covariates on the training folds.
        X_train_raw = [feat_train, cov_train];
        X_test_raw  = [feat_test,  cov_test];
        mu_X    = mean(X_train_raw, 1);
        sigma_X = std(X_train_raw, 0, 1);
        sigma_X(sigma_X == 0) = 1;
        X_train_std = (X_train_raw - mu_X) ./ sigma_X;
        X_test_std  = (X_test_raw  - mu_X) ./ sigma_X;

        % Class-balance weights inversely proportional to class size.
        n_pos = sum(y_train == 1);
        n_neg = sum(y_train == 0);
        if n_pos == 0 || n_neg == 0
            continue;
        end
        weights = ones(size(y_train));
        weights(y_train == 0) = n_pos / n_neg;

        try
            b = glmfit(X_train_std, y_train, 'binomial', 'weights', weights, 'constant', 'on');
        catch
            continue;
        end

        logit = [ones(size(X_test_std, 1), 1), X_test_std] * b;
        prob  = 1 ./ (1 + exp(-logit));

        if numel(unique(y_test)) < 2
            continue;
        end
        [~, ~, ~, auc_folds(fold)] = perfcurve(y_test, prob, 1);

        % Accuracy / sensitivity / specificity at the fixed probability cutoff.
        pred_class = prob >= cfg.class_threshold;
        TP = sum(pred_class == 1 & y_test == 1);
        FN = sum(pred_class == 0 & y_test == 1);
        FP = sum(pred_class == 1 & y_test == 0);
        TN = sum(pred_class == 0 & y_test == 0);

        acc_folds(fold) = (TP + TN) / numel(y_test);
        if TP + FN > 0
            sen_folds(fold) = TP / (TP + FN);
        end
        if TN + FP > 0
            spec_folds(fold) = TN / (TN + FP);
        end
    end

    valid_folds = ~isnan(auc_folds);
    if nnz(valid_folds) >= cfg.min_valid_folds
        auc_per_iter(iter)   = mean(auc_folds(valid_folds));
        edges_per_iter(iter) = mean(edges_folds(valid_folds));
        acc_per_iter(iter)   = mean(acc_folds(valid_folds), 'omitnan');
        sen_per_iter(iter)   = mean(sen_folds(valid_folds), 'omitnan');
        spec_per_iter(iter)  = mean(spec_folds(valid_folds), 'omitnan');
    else
        warning('Iteration %d: only %d/%d valid folds; this iteration is excluded.', ...
            iter, nnz(valid_folds), n_folds);
    end

    if mod(iter, 100) == 0
        fprintf('Completed %d/%d iterations\n', iter, n_iterations);
    end
end

%% Summary over valid iterations
valid_iters = ~isnan(auc_per_iter);
if nnz(valid_iters) < 2
    error('Fewer than 2 valid iterations (%d/%d); cannot summarize results.', ...
        nnz(valid_iters), n_iterations);
end

results.auc_per_iter   = auc_per_iter;
results.acc_per_iter   = acc_per_iter;
results.sen_per_iter   = sen_per_iter;
results.spec_per_iter  = spec_per_iter;
results.edges_per_iter = edges_per_iter;
results.rand_inds_all  = rand_inds_all;

results.mean_auc   = mean(auc_per_iter(valid_iters));
results.std_auc    = std(auc_per_iter(valid_iters));
results.mean_acc   = mean(acc_per_iter(valid_iters));
results.std_acc    = std(acc_per_iter(valid_iters));
results.mean_sen   = mean(sen_per_iter(valid_iters));
results.std_sen    = std(sen_per_iter(valid_iters));
results.mean_spec  = mean(spec_per_iter(valid_iters));
results.std_spec   = std(spec_per_iter(valid_iters));
results.mean_edges = mean(edges_per_iter(valid_iters));
results.std_edges  = std(edges_per_iter(valid_iters));

results.edge_selected_count     = edge_selected_count;
results.pos_edge_selected_count = pos_edge_selected_count;
results.neg_edge_selected_count = neg_edge_selected_count;

results.valid_iterations = nnz(valid_iters);
end
