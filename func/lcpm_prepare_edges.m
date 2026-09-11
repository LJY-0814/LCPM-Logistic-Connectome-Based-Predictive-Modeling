function [x_edges, triu_mask, n_roi] = lcpm_prepare_edges(x)
%LCPM_PREPARE_EDGES Vectorize connectivity matrices and apply Fisher z.
% Input:
%   x - 3-D ROI x ROI x n_subjects connectivity matrices
% Outputs:
%   x_edges   - n_edges x n_subjects matrix of Fisher-z transformed
%               upper-triangular values (one column per subject)
%   triu_mask - logical mask of the retained upper triangle
%   n_roi     - number of ROIs

if ndims(x) ~= 3 || size(x, 1) ~= size(x, 2)
    error('x must be a 3-D square connectivity array (ROI x ROI x subjects).');
end
n_roi     = size(x, 1);
n_subjects = size(x, 3);

triu_mask = logical(triu(ones(n_roi), 1));
x_edges   = zeros(nnz(triu_mask), n_subjects);

for i = 1:n_subjects
    r = x(:, :, i);
    x_edges(:, i) = r(triu_mask);
end
x_edges = atanh(x_edges);
end
