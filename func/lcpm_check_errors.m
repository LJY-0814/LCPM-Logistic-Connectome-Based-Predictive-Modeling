function lcpm_check_errors(x, y, folds, covariates)
%LCPM_CHECK_ERRORS Checks that input data are in a format usable by LCPM.
% Adapted from cpm_check_errors.m of the connectome-based predictive
% modeling toolbox (github.com/esfinn/cpm_tutorial), with additional
% checks for binary labels and covariates. x may be a 3-D
% ROI x ROI x subject connectivity array or a 2-D edge x subject matrix.

% Check that x data are in the required format - 2D or 3D
if (~ismatrix(x)) && (ndims(x) ~= 3)
    error('Data should have two or three dimensions');
end

% Check that x data contain more than one element
if size(x, 1) == 1
    error('Single feature detected.');
end

if size(x, ndims(x)) ~= size(y, 1)
    error('There are NOT the same number of subjects in the data and behavior variable');
end

% Check to make sure there are at least ten subjects in the input data
if size(x, ndims(x)) < 10
    warning('The LCPM code requires >10 subjects to function properly; sound results likely require >>10.');
end

% Check to make sure you have more subjects than folds
if size(x, ndims(x)) < folds
    warning('You must have more subjects than folds in your cross validation. Please check the help documentation.');
end

% Check whether x is symmetric across first two dimensions
if ndims(x) == 3
    if size(x, 1) ~= size(x, 2)
        warning('Please make sure, if intended, that data is an NxN connectivity matrix');
    end
end

% Check for nodes with values of 0 (missing nodes within a subject)
row_sum            = squeeze(sum(abs(x), 2));
zero_node          = sum(row_sum == 0);
zero_node_subjects = sum(zero_node > 0);
if zero_node_subjects > 0
    warning('Data: %d subjects have missing nodes. Please check your data.', zero_node_subjects);
end

% Check for Inf or NaN
if ~isempty(find(isinf(x), 1))
    warning('You have Inf values in your matrices. Please check your data.');
end

if ~isempty(find(isnan(x), 1))
    warning('You have NaNs in your matrices. Please check your data.');
end

%% Additional checks specific to the LCPM classifier

% Labels must be a binary vector
if ~isnumeric(y) || ~isvector(y)
    error('y must be a numeric vector.');
end
y = y(:);
if numel(unique(y(~isnan(y)))) ~= 2
    error('y must contain exactly two classes (binary labels, e.g. 0/1).');
end

% Covariates, if provided, must have one row per subject and no NaN/Inf
if nargin >= 4 && ~isempty(covariates)
    if size(covariates, 1) ~= numel(y)
        error('covariates must have one row per subject (%d rows expected, %d found).', ...
            numel(y), size(covariates, 1));
    end
    if any(isnan(covariates(:))) || any(isinf(covariates(:)))
        error('covariates contain NaN or Inf values. Please check your data.');
    end
end
end
