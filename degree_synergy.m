% Input weight matrices must have presynaptic neurons represented as rows,
% and postsynaptic neurons as columns.

function [corr_sender, corr_receiver] = degree_synergy(PID, weight_matrix, triplet_list, threshold)
    % Check inputs.
    if nargin < 2
        error('At least 2 input arguments are required.')
    elseif ~ismatrix(PID) || ~ismatrix(weight_matrix)
        error('Input datatype must be a matrix.')
    elseif size(PID,2) ~= 7
        error('Input PID values must a matrix with 7 columns.')
    elseif size(weight_matrix,1) ~= size(weight_matrix,2)
        error('Input adjacency matrix must be square.')
    end
    % Threshold.
    if nargin == 4
        if ~isscalar(threshold)
            error('Threshold must be a scalar.')
        elseif (threshold>1) || (threshold<0)
            error('Threshold must be between 0 and 1.')
        else
            ordered_weights = sort(weight_matrix(:));
            ordered_weights(ordered_weights==0) = [];
            threshold_value = ordered_weights(floor(threshold*size(ordered_weights,1)));
            % Binarize adjacency matrix
            weight_matrix(weight_matrix<threshold_value) = 0;
        clear threshold threshold_value ordered_weights
        end
    end
    weight_matrix(weight_matrix~=0) = 1; % Binarize adjacency matrix.
    indegree = sum(weight_matrix,1)';
    outdegree = sum(weight_matrix,2);
    if nargin < 3 || isempty(triplet_list)
        % Initialize nx3 matrix of all targeted non-zero neuron triplets, where column one indicates the target neuron.
        neuron_list = 1:size(weight_matrix,1);
        target_1 = nchoosek(neuron_list,3);
        target_2 = circshift(target_1,1,2);
        target_3 = circshift(target_1,-1,2);
        triplet_list = [target_1; target_2; target_3];
        clear target_1 target_2 target_3 neuron_list
    else
        if ~ismatrix(triplet_list)
            error('List of neuron triplets must be a matrix.')
        elseif size(triplet_list,2) ~= 3
            error('List of neuron triplets must have 3 columns.')
        end
    end
    [~,iAll] = intersect(PID(:,1:3), triplet_list, 'rows');
    clear triplet_list
    synergies = PID(iAll,4); % Extract synergy values.
    corr_sender = corrcoef([outdegree(PID(iAll,2)); outdegree(PID(iAll,3))], [synergies; synergies]); corr_sender = corr_sender(1,2);
    corr_receiver = corrcoef(indegree(PID(iAll,1)), synergies); corr_receiver = corr_receiver(1,2);
end