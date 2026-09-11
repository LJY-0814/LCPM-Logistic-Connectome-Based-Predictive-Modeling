function cfg = config()
%CONFIG All file paths and analysis parameters for the LCPM pipeline.
%   cfg = config() returns a struct with every path and threshold used by
%   main.m and the functions under func/. The thresholds are kept as
%   separate fields so the pipeline generalizes beyond the settings used
%   in the manuscript; change them here instead of editing the code.
%
% The default values reproduce the analyses reported in the manuscript:
% feature-selection p threshold 0.15, classification probability cutoff
% 0.40, high-frequency edge cutoff 90%, 10-fold cross-validation repeated
% 1,000 times.

%% File paths
% Paths are resolved relative to the project root (the folder containing
% this config.m), so the project can be moved without editing the code.
% Only data_root may need to point to your local copy of the dataset.
project_root = fileparts(mfilename('fullpath'));
data_root    = fullfile(project_root, '..', 'REST-meta-MDD');

cfg.data_dir       = fullfile(data_root, 'ROISignals MDD (all 2)');
cfg.roi_table_file = fullfile(data_root, 'SRDN ROI.xlsx');
cfg.label_file     = fullfile(data_root, 'y-SI MDD (all 2).xlsx');
cfg.output_dir     = fullfile(project_root, 'output');

%% Dataset
cfg.n_subjects = 183;   % must match the number of ROISignals_*.mat files

% Column indices in the label spreadsheet (no header row):
cfg.label_column      = 2;          % binary outcome: 1 = suicide risk, 0 = no risk
cfg.covariate_columns = [4 5 6 7];  % [HAMD total, sex, age, education]

%% Model parameters (thresholds are separate on purpose)
cfg.p_threshold     = 0.15;  % feature selection: point-biserial correlation p-value cutoff
cfg.class_threshold = 0.40;  % classification: probability cutoff for predicted labels
cfg.freq_threshold  = 0.90;  % consensus pattern: fraction of iterations x folds defining high-frequency edges

cfg.n_iterations    = 1000;  % repetitions of the k-fold cross-validation
cfg.n_folds         = 10;    % folds in each cross-validation iteration
cfg.min_valid_folds = 3;     % minimum number of valid folds for an iteration to be retained

%% Permutation test
cfg.permutation_seed = 0;    % seed for label shuffling; one permutation is run per
                             % original repeated-CV split, so the number of
                             % permutations equals cfg.n_iterations

%% Virtual lesion analysis
cfg.parallel    = false;     % use parfor over networks (requires Parallel Computing Toolbox)
cfg.fdr_correct = true;      % Benjamini-Hochberg FDR correction across networks
cfg.fdr_q       = 0.05;
end
