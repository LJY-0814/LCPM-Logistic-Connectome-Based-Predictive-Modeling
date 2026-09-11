function results = lcpm_permutation(x, y, covariates, rand_inds_all, observed_mean_auc, cfg)
%LCPM_PERMUTATION Label-permutation test for the LCPM classifier.
%   Inputs:
%     x          - n_edges x n_subjects matrix of Fisher-z connectivity
%                  values (from lcpm_prepare_edges)
%     y          - n_subjects x 1 binary labels (0/1)
%     covariates - n_subjects x n_covariates matrix
%     rand_inds_all - n_subjects x cfg.n_iterations matrix of predefined
%                     splits from lcpm_repeated_cv
%     observed_mean_auc - mean AUC from the real repeated-CV run
%     cfg        - parameter struct from config.m
%
%   Output:
%     results - struct with permutation AUCs, p value, and diagnostics

n_perm = cfg.n_iterations;
n_subjects = size(x, 2);
n_folds = cfg.n_folds;

if size(rand_inds_all, 2) < n_perm
    error('rand_inds_all has %d columns but cfg.n_iterations = %d; need one split per permutation.', ...
        size(rand_inds_all, 2), n_perm);
end

base_size = floor(n_subjects / n_folds);
remainder = n_subjects - base_size * n_folds;

perm_aucs = nan(n_perm, 1);
fprintf('Running %d permutations...\n', n_perm);
rng(cfg.permutation_seed, 'twister');

for i = 1:n_perm
    y_shuffled = y(randperm(n_subjects));
    randinds = rand_inds_all(:, i);

    auc_folds = nan(n_folds, 1);
    start_idx = 1;

    for fold = 1:n_folds
        test_size = base_size + (fold <= remainder);
        test_pos  = start_idx : (start_idx + test_size - 1);
        start_idx = start_idx + test_size;
        train_pos = [1:test_pos(1) - 1, test_pos(end) + 1:n_subjects];

        test_inds  = randinds(test_pos);
        train_inds = randinds(train_pos);

        X_train = x(:, train_inds).';
        y_train = y_shuffled(train_inds);
        cov_train = covariates(train_inds, :);
        X_test  = x(:, test_inds).';
        y_test  = y_shuffled(test_inds);
        cov_test = covariates(test_inds, :);

        % Feature selection on the training folds only.
        [rho, pval] = corr(X_train, y_train, 'Type', 'Pearson');
        pos_idx = find(pval < cfg.p_threshold & rho > 0);
        neg_idx = find(pval < cfg.p_threshold & rho < 0);

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
    end

    valid_folds = ~isnan(auc_folds);
    if nnz(valid_folds) >= cfg.min_valid_folds
        perm_aucs(i) = mean(auc_folds(valid_folds));
    else
        warning('Permutation %d had only %d/%d valid folds; recorded as NaN.', ...
            i, nnz(valid_folds), n_folds);
    end

    if mod(i, 100) == 0
        fprintf('Permutation progress: %d/%d\n', i, n_perm);
    end
end

valid_perms = ~isnan(perm_aucs);
perm_p = sum(perm_aucs(valid_perms) >= observed_mean_auc) / nnz(valid_perms);

fprintf('\nPermutation test (%d permutations):\n', n_perm);
fprintf('Observed mean AUC        = %.4f\n', observed_mean_auc);
fprintf('Permutation AUC          = %.4f +/- %.4f\n', mean(perm_aucs(valid_perms)), std(perm_aucs(valid_perms)));
fprintf('One-tailed p value       = %.4f\n', perm_p);
if nnz(~valid_perms) > 0
    warning('%d/%d permutations had no valid folds and were excluded from the p value.', ...
        nnz(~valid_perms), n_perm);
end

results.perm_aucs      = perm_aucs;
results.perm_p         = perm_p;
results.observed_auc   = observed_mean_auc;
results.n_permutations = n_perm;
results.n_valid        = nnz(valid_perms);
end
