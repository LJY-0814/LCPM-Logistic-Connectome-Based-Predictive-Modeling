%% main.m
% Entry point of the LCPM pipeline.
%
% Steps:
%   1. Load connectivity data, labels and covariates   (load_data)
%   2. Validate inputs                               (lcpm_check_errors)
%   3. Vectorize matrices and apply Fisher z         (lcpm_prepare_edges)
%   4. Repeated k-fold cross-validation              (lcpm_repeated_cv)
%   5. Label-permutation test on the same splits     (lcpm_permutation)
%   6. Consensus-pattern characterization            (lcpm_characterize_pattern)
%   7. Virtual lesion (leave-one-network-out)        (lcpm_virtual_lesion)
%
% All plotting code has been removed; every result is written as a table
% or a MAT file under cfg.output_dir. Edit config.m before running.

clear; clc;

cfg = config();
if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

%% 1. Data loading
[data, roi_info] = load_data(cfg);

%% 2. Input validation
lcpm_check_errors(data.x, data.y, cfg.n_folds, data.covariates);

%% 3. Fisher-z vectorization
[x_edges, ~, n_roi] = lcpm_prepare_edges(data.x);
fprintf('Loaded %d subjects, %d ROIs, %d edges per subject.\n', ...
    cfg.n_subjects, n_roi, size(x_edges, 1));

%% 4. Repeated k-fold cross-validation
fprintf('\n===== Repeated %d-fold cross-validation (%d iterations, p = %.3f, classification threshold = %.2f) =====\n', ...
    cfg.n_folds, cfg.n_iterations, cfg.p_threshold, cfg.class_threshold);
cv_results = lcpm_repeated_cv(x_edges, data.y, data.covariates, cfg);

fprintf('\n===== Cross-validation results =====\n');
fprintf('Mean AUC          = %.4f +/- %.4f\n', cv_results.mean_auc, cv_results.std_auc);
fprintf('Mean accuracy     = %.4f +/- %.4f\n', cv_results.mean_acc, cv_results.std_acc);
fprintf('Mean sensitivity  = %.4f +/- %.4f\n', cv_results.mean_sen, cv_results.std_sen);
fprintf('Mean specificity  = %.4f +/- %.4f\n', cv_results.mean_spec, cv_results.std_spec);
fprintf('Mean edge count   = %.2f +/- %.2f\n', cv_results.mean_edges, cv_results.std_edges);
fprintf('Valid iterations  = %d / %d\n', cv_results.valid_iterations, cfg.n_iterations);

summary = table(cfg.p_threshold, cfg.class_threshold, ...
    cv_results.mean_auc, cv_results.std_auc, ...
    cv_results.mean_acc, cv_results.std_acc, ...
    cv_results.mean_sen, cv_results.std_sen, ...
    cv_results.mean_spec, cv_results.std_spec, ...
    cv_results.mean_edges, cv_results.std_edges, ...
    cv_results.valid_iterations, ...
    'VariableNames', {'p_threshold', 'classification_threshold', ...
    'mean_AUC', 'std_AUC', 'mean_accuracy', 'std_accuracy', ...
    'mean_sensitivity', 'std_sensitivity', 'mean_specificity', 'std_specificity', ...
    'mean_edges', 'std_edges', 'valid_iterations'});
writetable(summary, fullfile(cfg.output_dir, 'cv_summary.xlsx'));
fprintf('Summary saved to %s\n', fullfile(cfg.output_dir, 'cv_summary.xlsx'));

%% 5. Permutation test (one permutation per original split, shuffled labels)
perm_results = lcpm_permutation(x_edges, data.y, data.covariates, ...
    cv_results.rand_inds_all, cv_results.mean_auc, cfg);
perm_table = table(perm_results.observed_auc, mean(perm_results.perm_aucs, 'omitnan'), ...
    std(perm_results.perm_aucs, 'omitnan'), perm_results.perm_p, perm_results.n_valid, ...
    'VariableNames', {'observed_mean_AUC', 'perm_mean_AUC', 'perm_std_AUC', ...
    'one_tailed_p', 'n_valid_permutations'});
writetable(perm_table, fullfile(cfg.output_dir, 'permutation_results.xlsx'));

%% 6. Consensus-pattern characterization (high-frequency edges, node degree,network-level distribution)
pattern_results = lcpm_characterize_pattern(cv_results, roi_info, cfg);

%% 7. Virtual lesion analysis (same splits, one network removed at a time)
lesion_results = lcpm_virtual_lesion(x_edges, data.y, data.covariates, ...
    roi_info.network_names, cv_results.rand_inds_all, cv_results.auc_per_iter, cfg);

%% Save everything for downstream use
save(fullfile(cfg.output_dir, 'lcpm_results.mat'), ...
    'cfg', 'data', 'roi_info', 'cv_results', 'perm_results', 'pattern_results', 'lesion_results');
fprintf('\nAll results saved under %s\n', cfg.output_dir);
