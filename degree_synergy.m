% Given a weight matrix and a list of PID values, return correlation
% between source neuron out-degree and synergy, as well as target in-degree
% and synergy. A list of all neuron in and out degrees may also be
% returned.
%
% Optionally, input a list of neuron entropies over which to normalize
% synergy values; input a list of neuron triplets over which to calculate
% correlation; and a threshold between 0 and 1. E.g. if a threshold of 0.25
% is given, set the bottom 25% of weights in the weight matrix to zero.
%
% Input weight matrices must have presynaptic neurons represented as rows,
% and postsynaptic neurons as columns.

function [corr_source_out, corr_target_in, indegree, outdegree] = degree_synergy(PID, weight_matrix, entropy, triplet_list, threshold)

    %% Check inputs.
    if nargin < 2
        error('At least 2 input arguments are required: PID values and a weight matrix.')
    elseif ~ismatrix(PID) || ~ismatrix(weight_matrix)
        error('Input datatype must be a matrix.')
    elseif size(PID,2) ~= 7
        error('Input PID values must a matrix with 7 columns, e.g. output of TE_PID.m.')
    elseif size(weight_matrix,1) ~= size(weight_matrix,2)
        error('Input adjacency matrix must be square.')
    end
    nNeuron = size(weight_matrix,1);
    
    %% Optionally, threshold weight matrix.
    if nargin == 5
        if ~isscalar(threshold)
            error('Threshold must be a scalar.')
        elseif (threshold>1) || (threshold<0)
            error('Threshold must be between 0 and 1.')
        else
            ordered_weights = sort(weight_matrix(:));
            ordered_weights(ordered_weights==0) = [];
            threshold_value = ordered_weights(floor(threshold*size(ordered_weights,1)));
            weight_matrix(weight_matrix<threshold_value) = 0;
        clear threshold threshold_value ordered_weights
        end
    end
    
    %% Create binary adjacency matrix from weight matrix, then calculate in and out degrees.
    weight_matrix(weight_matrix~=0) = 1; % Binarize adjacency matrix.
    indegree = sum(weight_matrix,1)'; % Returns a 1-dimensional list of neuron in-degrees, the indices of which correspond to indices of the weight matrix.
    outdegree = sum(weight_matrix,2); % Similarly, returns neuron out-degrees.
    
    %% Optionally, specify a subset of neuron triplets over which to calculate synergy-degree correlation.
    if nargin < 4 || isempty(triplet_list)
        iAll = true(size(PID,1),1);
    else
        % Check given input triplet list.
        if ~ismatrix(triplet_list)
            error('List of neuron triplets must be a matrix.')
        elseif size(triplet_list,2) ~= 3
            error('List of neuron triplets must have 3 columns, specifying target neuron index, and two source neuron indices.')
        elseif any(unique(triplet_list) > size(input_data,2))
            error('Neuron indices in given triplet list must not be greater than total number of neurons in input dataset.')
        else
            [~,iAll] = intersect(PID(:,1:3), triplet_list, 'rows');
        end
    end
    
    %% Optionally, normalize synergy values by entropy of target neuron.
    % Only triplets present in both input PID values and triplet list will
    % be considered. Zero entropy cases are thus ignored if input PID
    % values are outputs of TE_PID.m.
    if (nargin > 2) && ~isempty(entropy)
        if size(entropy,2) ~= 2
            error('Input entropies must have 2 columns, denoting neuron index and entropy.')
        % Ensure length of entropy list is equal to number of neurons.
        elseif size(entropy,1) ~= nNeuron
            extended_entropy = [(1:nNeuron)' zeros(nNeuron,1)];
            for i = 1:nNeuron
                if sum(entropy(:,1)==i) ~= 0
                    extended_entropy(i,2) = entropy(entropy(:,1)==i,2);
                end
            end
            entropy = extended_entropy; % Rename.
            clear extended_entropy
        end
        entropies = entropy(PID(iAll,1),2);
        if sum(entropies==0) ~= 0
            iNonzeroEntropy = entropy(PID(:,1),2)~=0;
            iAll = iAll & iNonzeroEntropy;
        end
        synergies = PID(iAll,4);
        entropies = entropy(PID(iAll,1),2);
        synergies = synergies ./ entropies;
    else
        synergies = PID(iAll,4);
    end
    
    %% Return correlations between synergy and target in-degree, as well as synergy and source out-degree.
    corr_source_out = corrcoef([outdegree(PID(iAll,2)); outdegree(PID(iAll,3))], [synergies; synergies]); corr_source_out = corr_source_out(2);
    corr_target_in = corrcoef(indegree(PID(iAll,1)), synergies); corr_target_in = corr_target_in(2);
    
end