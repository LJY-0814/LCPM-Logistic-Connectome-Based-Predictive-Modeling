function [data, roi_info] = load_data(cfg)
%LOAD_DATA Load connectivity matrices, labels and covariates.
% This loader is specific to the REST-meta-MDD ROISignals files used in the
% manuscript (one ROISignals_S*.mat file per subject, ROI column indices in
% an Excel table). To apply LCPM to another dataset, replace this function
% with your own; the downstream pipeline only needs:
%   data.x           - ROI x ROI x n_subjects connectivity matrices (raw
%                      Pearson correlations)
%   data.y           - n_subjects x 1 binary labels (0/1)
%   data.covariates  - n_subjects x n_covariates matrix
%   roi_info         - struct with fields roi_numbers, region_names,
%                      network_names (one entry per ROI)
%
%   cfg.roi_table_file - ROI table with at least 8 columns; only three are
%       used: column 1 = ROI column indices (which ROISignals columns to
%       keep), column 2 = region names, column 8 = network names (lowercased
%       and trimmed). Columns 3-7 are ignored.

%% ROI definition
roi_table = readtable(cfg.roi_table_file);
if width(roi_table) < 8
    error(['ROI table %s must have at least 8 columns: ', ...
           'column 1 = ROI column indices, column 2 = region names, ', ...
           'column 8 = network names.'], cfg.roi_table_file);
end
roi_columns   = roi_table{:, 1};
region_names  = string(roi_table{:, 2});
network_names = lower(strtrim(string(roi_table{:, 8})));
n_roi         = numel(roi_columns);
fprintf('Loaded %d ROIs from %s\n', n_roi, cfg.roi_table_file);

%% Subject files
all_mat_files = dir(fullfile(cfg.data_dir, '*.mat'));
if isempty(all_mat_files)
    error('No .mat files found in %s', cfg.data_dir);
end

% Keep only files matching the expected naming scheme ROISignals_S*.mat and
% record the numeric tokens so subjects can be sorted deterministically.
file_nums = nan(numel(all_mat_files), 3);
keep      = false(numel(all_mat_files), 1);
for i = 1:numel(all_mat_files)
    tokens = regexp(all_mat_files(i).name, 'ROISignals_S(\d+)-(\d+)-(\d+)\.mat', 'tokens');
    if ~isempty(tokens)
        file_nums(i, :) = [str2double(tokens{1}{1}), str2double(tokens{1}{2}), str2double(tokens{1}{3})];
        keep(i)         = true;
    end
end
all_mat_files = all_mat_files(keep);
file_nums     = file_nums(keep, :);
if numel(all_mat_files) ~= cfg.n_subjects
    error('Expected %d subject files but found %d matching ROISignals_S*.mat in %s.', ...
        cfg.n_subjects, numel(all_mat_files), cfg.data_dir);
end

% Sort subjects by site id and subject id.
[~, sort_idx] = sortrows(file_nums, [1, 3]);
all_mat_files = all_mat_files(sort_idx);

%% Connectivity matrices
x = zeros(n_roi, n_roi, cfg.n_subjects);
for i = 1:cfg.n_subjects
    fprintf('Processing file %d/%d: %s\n', i, cfg.n_subjects, all_mat_files(i).name);
    s = load(fullfile(cfg.data_dir, all_mat_files(i).name), 'ROISignals');
    if ~isfield(s, 'ROISignals')
        error('File %s does not contain a ROISignals variable.', all_mat_files(i).name);
    end
    conn = corrcoef(s.ROISignals(:, roi_columns));
    if size(conn, 1) ~= n_roi
        error('File %s produced a %d x %d matrix, expected %d x %d.', ...
            all_mat_files(i).name, size(conn, 1), size(conn, 2), n_roi, n_roi);
    end
    x(:, :, i) = conn;
end

%% Labels and covariates
label_table = readtable(cfg.label_file, 'ReadVariableNames', false);
if size(label_table, 1) ~= cfg.n_subjects
    error('Label file %s has %d rows but %d subjects were loaded.', ...
        cfg.label_file, size(label_table, 1), cfg.n_subjects);
end
needed_cols = [cfg.label_column, cfg.covariate_columns];
if max(needed_cols) > width(label_table)
    error('Label file %s has %d columns; column index %d is required.', ...
        cfg.label_file, width(label_table), max(needed_cols));
end
y          = double(label_table{:, cfg.label_column});
covariates = label_table{:, cfg.covariate_columns};

fprintf('Total subjects: %d (label 1: %d, label 0: %d)\n', ...
    cfg.n_subjects, sum(y == 1), sum(y == 0));

data     = struct('x', x, 'y', y(:), 'covariates', covariates);
roi_info = struct('roi_numbers', roi_columns, ...
                  'region_names', region_names, ...
                  'network_names', network_names);
end
