% Given a functional matrix and a synaptic matrix, this function returns
% the Pearson correlation coefficient between elements with identical
% indices in the two matrices. Only non-zero functional connections (and
% their corresponding synaptic connections) are considered.

function [zscores, corr_func_synp] =  compare_matrix(functional_matrix, synaptic_matrix, threshold)
    if ~ismatrix(functional_matrix) || ~ismatrix(synaptic_matrix)
        error('Input weight matrices must be matrices.')
    elseif (size(functional_matrix,1)~=size(synaptic_matrix,1)) || (size(functional_matrix,2)~=size(synaptic_matrix,2))
        error('Input weight matrices must have the same dimensions.')
    elseif size(functional_matrix,1)~=size(functional_matrix,2)
        error('Input weight matrices must be square.')
    end
    func_weights = functional_matrix(:);
    synp_weights = synaptic_matrix(:);
    synp_weights(func_weights==0) = [];
    func_weights(func_weights==0) = [];
    func_zscore = zscore(func_weights);
    synp_zscore = zscore(synp_weights);
    if nargin == 3
        if ~isscalar(threshold)
            error('Threshold must be a scalar.')
        elseif (threshold>1) || (threshold<0)
            error('Threshold must be between 0 and 1.')
        else
            ordered_weights = sort(func_weights);
            threshold_value = ordered_weights(floor(threshold*size(ordered_weights,1)));
            synp_zscore(func_weights<threshold_value) = [];
            func_zscore(func_weights<threshold_value) = [];
        end
    end
    zscores = [func_zscore synp_zscore];
    corr_func_synp = corrcoef(func_zscore, synp_zscore); corr_func_synp = corr_func_synp(2);
end