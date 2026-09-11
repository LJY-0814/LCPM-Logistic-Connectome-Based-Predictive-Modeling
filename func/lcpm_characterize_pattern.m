function results = lcpm_characterize_pattern(selection, roi_info, cfg)
%LCPM_CHARACTERIZE_PATTERN Consensus-pattern characterization.
% Characterizes the predictive pattern that is stable across the repeated
% cross-validation runs:
%   1. High-frequency edges: edges selected in at least
%      cfg.freq_threshold of all iterations x folds form the consensus
%      network.
%   2. Node degree: number of high-frequency edges incident to each ROI.
%   3. Network-level distribution: number of high-frequency edges
%      within and between functional networks.
%
% Inputs:
%   selection  - struct from lcpm_repeated_cv with fields
%                edge_selected_count / pos_edge_selected_count /
%                neg_edge_selected_count
%   roi_info   - struct from load_data with fields roi_numbers,
%                region_names, network_names (one entry per ROI)
%   cfg        - parameter struct from config.m
%
% Tables are written to cfg.output_dir.

%% Validate inputs
n_roi = numel(roi_info.region_names);
if numel(roi_info.network_names) ~= n_roi || numel(roi_info.roi_numbers) ~= n_roi
    error('roi_info fields must have one entry per ROI (%d).', n_roi);
end
expected_edges = n_roi * (n_roi - 1) / 2;
for f = {'edge_selected_count', 'pos_edge_selected_count', 'neg_edge_selected_count'}
    if ~isfield(selection, f{1}) || numel(selection.(f{1})) ~= expected_edges
        error('selection.%s must have one entry per edge (%d).', f{1}, expected_edges);
    end
end

triu_mask = triu(true(n_roi), 1);
[roi_i, roi_j] = ind2sub([n_roi, n_roi], find(triu_mask));

%% 1. High-frequency edges
total_folds = cfg.n_iterations * cfg.n_folds;
cutoff      = cfg.freq_threshold * total_folds;

selected_edges = find(selection.edge_selected_count >= cutoff);
n_selected     = numel(selected_edges);
n_pos_sel      = nnz(selection.pos_edge_selected_count(selected_edges) >= cutoff);
n_neg_sel      = nnz(selection.neg_edge_selected_count(selected_edges) >= cutoff);

fprintf('High-frequency edges (>= %.0f%% of %d folds): %d\n', ...
    100 * cfg.freq_threshold, total_folds, n_selected);
fprintf('  selected as positive in >= %.0f%% of folds: %d\n', 100 * cfg.freq_threshold, n_pos_sel);
fprintf('  selected as negative in >= %.0f%% of folds: %d\n', 100 * cfg.freq_threshold, n_neg_sel);

region_names = string(roi_info.region_names);

edge_info = table();
edge_info.EdgeIndex      = selected_edges;
edge_info.ROI1           = roi_i(selected_edges);
edge_info.ROI2           = roi_j(selected_edges);
edge_info.ROINumber1     = roi_info.roi_numbers(edge_info.ROI1);
edge_info.ROINumber2     = roi_info.roi_numbers(edge_info.ROI2);
edge_info.Region1        = region_names(edge_info.ROI1);
edge_info.Region2        = region_names(edge_info.ROI2);
edge_info.SelectedCount  = selection.edge_selected_count(selected_edges);
edge_info.Frequency      = selection.edge_selected_count(selected_edges) / total_folds;

% Dominant sign of each high-frequency edge (counts across all folds).
edge_info.Sign = repmat("", n_selected, 1);
for e = 1:n_selected
    idx = selected_edges(e);
    if selection.pos_edge_selected_count(idx) >= selection.neg_edge_selected_count(idx)
        edge_info.Sign(e) = "positive";
    else
        edge_info.Sign(e) = "negative";
    end
end

out_file = fullfile(cfg.output_dir, 'selected_edges_high_frequency.xlsx');
writetable(edge_info, out_file);
fprintf('High-frequency edge table saved to %s\n', out_file);

%% 2. Node degree
node_ids = [roi_i(selected_edges); roi_j(selected_edges)];
[u_nodes, ~, node_grp] = unique(node_ids);
degree = accumarray(node_grp, 1);

degree_table = table();
degree_table.ROIIndex  = u_nodes;
degree_table.ROINumber = roi_info.roi_numbers(u_nodes);
degree_table.Region    = region_names(u_nodes);
degree_table.Degree    = degree;
degree_table           = sortrows(degree_table, 'Degree', 'descend');

fprintf('Top 10 regions by degree in the consensus network:\n');
disp(degree_table(1:min(10, height(degree_table)), :));

out_file = fullfile(cfg.output_dir, 'node_degree_high_frequency.xlsx');
writetable(degree_table, out_file);
fprintf('Node degree table saved to %s\n', out_file);

%% 3. Network-level distribution
network_names = lower(strtrim(string(roi_info.network_names)));
network_names(ismissing(network_names) | network_names == "") = "__missing__";
networks = unique(network_names);
networks = networks(networks ~= "__missing__");
n_nets   = numel(networks);

net_i = network_names(roi_i(selected_edges));
net_j = network_names(roi_j(selected_edges));
edge_net_pairs = sort([net_i, net_j], 2); 

net_count = zeros(n_nets, n_nets);

for e = 1:n_selected
    idx1 = find(networks == edge_net_pairs(e, 1));
    idx2 = find(networks == edge_net_pairs(e, 2));
    if isempty(idx1) || isempty(idx2)
        continue;
    end
    net_count(idx1, idx2) = net_count(idx1, idx2) + 1;
    if idx1 ~= idx2
        net_count(idx2, idx1) = net_count(idx2, idx1) + 1;
    end
end

% Inter- vs intra-network split of the consensus edges 
n_intra = sum(diag(net_count));
if n_selected > 0
    fprintf('Intra-network high-frequency edges: %d (%.1f%%)\n', ...
        n_intra, 100 * n_intra / n_selected);
end

net_names_cell = cellstr(networks);
fprintf('High-frequency edge counts between networks:\n');
disp(array2table(net_count, 'VariableNames', net_names_cell, 'RowNames', net_names_cell));

out_file = fullfile(cfg.output_dir, 'network_edge_count.xlsx');
writetable(array2table(net_count, 'VariableNames', net_names_cell, 'RowNames', net_names_cell), ...
    out_file, 'WriteRowNames', true);

%% Collect outputs
results.edge_info    = edge_info;
results.degree_table = degree_table;
results.networks     = networks;
results.net_count    = net_count;
results.n_selected   = n_selected;
results.total_folds  = total_folds;
end
