% This function finds the recruitment matrix, i.e. the intersection of the
% functional matrix and the synaptic connectivity matrix.

function [recruitment_triplets, recruitment_matrix] = recruitment(functional_matrix, synaptic_matrix, threshold)
    if ~ismatrix(functional_matrix) || ~ismatrix(synaptic_matrix)
        error('Input weight matrices must be matrices.')
    elseif size(functional_matrix)~=size(synaptic_matrix)
        error('Input weight matrices must have the same dimensions.')
    elseif size(functional_matrix,1)~=size(functional_matrix,2)
        error('Input weight matrices must be square.')
    end
    if nargin == 3
        if ~isscalar(threshold)
            error('Threshold must be a scalar.')
        elseif (threshold>1) || (threshold<0)
            error('Threshold must be between 0 and 1.')
        else
            ordered_weights = synaptic_matrix(:);
            ordered_weights(ordered_weights==0) = [];
            ordered_weights = sort(ordered_weights);
            threshold_value = ordered_weights(floor(threshold*size(ordered_weights,1)));
            synaptic_matrix(synaptic_matrix<threshold_value) = 0;
        end
    end
    % Initialize recruitment network weight matrix.
    recruitment_matrix = functional_matrix;
    % Intersect functional matrix with synaptic matrix.
    recruitment_matrix(synaptic_matrix==0) = 0;
    % Initialize nx3 matrix of all targeted non-zero neuron triplets.
    length_vector = 1:size(functional_matrix,2);
    target_1 = nchoosek(length_vector,3);
    target_2 = circshift(target_1,1,2);
    target_3 = circshift(target_1,-1,2);
    all_triplets_list = [target_1; target_2; target_3];
    clear target_1
    clear target_2
    clear target_3
    % Initialize list of targeted neuron triplets of the recruitment network.
    syms n
    recruitment_triplets = zeros(double(symsum(nchoosek(n,2), n, 2, size(length_vector,2)-1)), 3);
    clear length_vector
    % Initialize dummy variable to indicate rows of output matrix.
    row_index = 1;
    for i = all_triplets_list'
        if (recruitment_matrix(i(1), i(2))>0) && (recruitment_matrix(i(1), i(3))>0)
            recruitment_triplets(row_index,:) = i'; % Write to row of output matrix indicated by row_index.
            row_index = row_index+1; % Increment dummy variable by 1.
        end
    end
    recruitment_triplets(recruitment_triplets(:,1)==0,:) = []; % Remove zeroes.
end