function results = lcpm_virtual_lesion(x, y, covariates, network_names, rand_inds_all, whole_aucs, cfg)
%LCPM_VIRTUAL_LESION Leave-one-network-out (virtual lesion) analysis.
%   results = lcpm_virtual_lesion(x, y, covariates, network_names, ...
%                                 rand_inds_all, whole_aucs, cfg)
%
% For each functional network, all edges involving at least one ROI of that
% network are removed and the full LCPM pipeline is re-run on the exact
% same data splits as the whole-brain analysis (rand_inds_all), so any
% performance change is attributable solely to the removed edges.
%
% the result table is written tovirtual_lesion_results.xlsx in cfg.output_dir.

%% Validate inputs
network_names = string(network_names(:).');
network_names(ismissing(network_names) | network_names == "") = "__missing__";
networks = unique(network_names);
networks = networks(networks ~= "__missing__");
n_net = numel(networks);
if n_net == 0
    error('No valid network labels found.');
end

n_roi = numel(network_names);
expected_edges = n_roi * (n_roi - 1) / 2;
if size(x, 1) ~= expected_edges
    error('x has %d edges but %d ROIs imply %d upper-triangular edges.', ...
        size(x, 1), n_roi, expected_edges);
end
if size(x, 2) ~= numel(y)
    error('x must have one column per subject.');
end

triu_mask = triu(true(n_roi), 1);
[roi_i, roi_j] = ind2sub([n_roi, n_roi], find(triu_mask));

%% Whole-brain reference AUCs on the same splits
if isempty(whole_aucs)
    fprintf('No whole-brain AUCs supplied; evaluating the whole-brain model on the same splits...\n');
    whole = lcpm_repeated_cv(x, y, covariates, cfg, 'RandInds', rand_inds_all);
    whole_aucs = whole.auc_per_iter;
end
if numel(whole_aucs) ~= cfg.n_iterations
    error('whole_aucs must have one entry per iteration (%d).', cfg.n_iterations);
end

%% Leave-one-network-out runs
net_aucs = nan(cfg.n_iterations, n_net);

if cfg.parallel
    parfor n = 1:n_net
        net_aucs(:, n) = lesioned_run(x, y, covariates, network_names, ...
            networks(n), roi_i, roi_j, rand_inds_all, cfg);
    end
else
    for n = 1:n_net
        fprintf('Virtual lesion: network %d/%d (%s)\n', n, n_net, networks(n));
        net_aucs(:, n) = lesioned_run(x, y, covariates, network_names, ...
            networks(n), roi_i, roi_j, rand_inds_all, cfg);
    end
end

%% Paired Wilcoxon signed-rank tests against the whole-brain model
valid_iter     = ~isnan(whole_aucs);
valid_whole    = whole_aucs(valid_iter);
mean_whole_auc = mean(valid_whole);

mean_net_auc = nan(1, n_net);
std_net_auc  = nan(1, n_net);
p_vs_whole   = nan(n_net, 1);
significantly_worse = false(n_net, 1);

for n = 1:n_net
    data = net_aucs(:, n);
    pair = valid_iter & ~isnan(data);
    mean_net_auc(n) = mean(data(pair), 'omitnan');
    std_net_auc(n)  = std(data(pair), [], 'omitnan');
    if nnz(pair) < 5
        continue;
    end
    % One-tailed: is the whole-brain AUC greater than the lesioned AUC?
    p_vs_whole(n) = signrank(valid_whole(pair(valid_iter)), data(pair), 'tail', 'right');
    significantly_worse(n) = p_vs_whole(n) < 0.05 && mean_net_auc(n) < mean_whole_auc;
end

fdr_significant = false(n_net, 1);
if cfg.fdr_correct
    fdr_significant = benjamini_hochberg(p_vs_whole, cfg.fdr_q);
end

%% Report and save
T = table(networks', mean_net_auc', std_net_auc', p_vs_whole, significantly_worse, ...
    'VariableNames', {'Network', 'Mean_AUC', 'Std_AUC', 'p_Wilcoxon', 'Significantly_Worse'});
if cfg.fdr_correct
    T.FDR_Significant = fdr_significant;
end
disp(T);

out_file = fullfile(cfg.output_dir, 'virtual_lesion_results.xlsx');
writetable(T, out_file);
fprintf('Virtual lesion results saved to %s\n', out_file);

results.networks            = networks;
results.net_aucs            = net_aucs;
results.whole_aucs          = whole_aucs;
results.mean_net_auc        = mean_net_auc;
results.std_net_auc         = std_net_auc;
results.mean_whole_auc      = mean_whole_auc;
results.p_vs_whole          = p_vs_whole;
results.significantly_worse = significantly_worse;
results.fdr_significant     = fdr_significant;
results.table               = T;
end

function aucs = lesioned_run(x, y, covariates, network_names, network, roi_i, roi_j, rand_inds_all, cfg)
% Re-run the LCPM with all edges touching ROIs of the given network removed.
roi_in_net = network_names == network;
keep       = ~(roi_in_net(roi_i) | roi_in_net(roi_j));
if ~any(keep)
    warning('Removing network %s removes all edges; AUC set to NaN.', network);
    aucs = nan(cfg.n_iterations, 1);
    return;
end
r    = lcpm_repeated_cv(x, y, covariates, cfg, 'RandInds', rand_inds_all, 'EdgeMask', keep);
aucs = r.auc_per_iter;
end

function h = benjamini_hochberg(pvals, q)
% Benjamini-Hochberg FDR correction.
h      = false(size(pvals));
ok     = ~isnan(pvals);
p      = pvals(ok);
m      = numel(p);
[~, order] = sort(p);
thresh = (1:m)' * q / m;
k      = find(p(order) <= thresh, 1, 'last');
if ~isempty(k)
    h_sorted           = false(m, 1);
    h_sorted(order(1:k)) = true;
    h(ok)              = h_sorted;
end
end
