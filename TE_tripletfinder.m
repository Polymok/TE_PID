% Given a matrix whose columns contain each time-series and a time-delay,
% this function returns a 3-column matrix whose columns indicate neuron
% indexes in increasing order:
% target | source1 | source2 
% satisfying TE(source1->target)>TE(target->source1) and
% TE(source2->target)>TE(target->source2) at the given time-delay.
%
% This function is designed to filter neuron triplets for transfer entropy
% partial information decomposition calculation.

function [functional_triplets, functional_matrix, threshold_matrix] = TE_tripletfinder(input_timeseries, delay, threshold)
    % Check if input formats are acceptable.
    if ~ismatrix(input_timeseries)
        error('Input time-series must be a matrix.')
    elseif ~isscalar(delay)
        error('Input time-delay must be a scalar.')
    elseif delay<1
        error('Input time-delay must be at least 1.')
    end
    % Check if a single time-series is contained in a column.
    if size(input_timeseries,1) < size(input_timeseries,2)
        str = input('Input matrix has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            input_timeseries = input_timeseries';
        end
    end
    % Initialize weight matrix of zeroes whose rows are source neuron
    % indices and columns are target neuron indices. For a given neuron
    % pair ij and a time-delay, if the greater transfer entropy obtains
    % when i lags behind j, we record TE(i->j) in element (i,j) of the
    % directional matrix. Conversely, if greater transfer entropy obtains
    % when j lags behind i, we record TE(j->i in element (j,i) of the
    % directional matrix.
    functional_matrix = zeros(size(input_timeseries,2), size(input_timeseries,2));
    % Initialize list of all undirected neuron pairs.
    length_vector = 1:size(input_timeseries,2);
    neuron_pairs = nchoosek(length_vector,2);
    for i = 1:size(neuron_pairs,1)
        [i_to_j, ~] = TE(input_timeseries(:,neuron_pairs(i,2)), input_timeseries(:,neuron_pairs(i,1)), delay);
        [j_to_i, ~] = TE(input_timeseries(:,neuron_pairs(i,1)), input_timeseries(:,neuron_pairs(i,2)), delay);
        if i_to_j > j_to_i
            functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) = i_to_j;
        elseif i_to_j < j_to_i
            functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1)) = j_to_i;
        elseif i_to_j==0 && j_to_i==0
        elseif i_to_j==j_to_i
            disp('Transfer entropy in both directions are equal.')
            functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) = i_to_j;
            functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1)) = i_to_j;
        end
    end
    clear neuron_pairs
    % Initialze weight matrix thresholded at given threshold.
    threshold_matrix = functional_matrix;
    if nargin == 3
        if ~isscalar(threshold)
            error('Threshold must be a scalar.')
        elseif (threshold>1) || (threshold<0)
            error('Threshold must be between 0 and 1.')
        else
            ordered_weights = functional_matrix(:);
            ordered_weights(ordered_weights==0) = [];
            ordered_weights = sort(ordered_weights);
            threshold_value = ordered_weights(floor(threshold*size(ordered_weights,1)));
            threshold_matrix(threshold_matrix<threshold_value) = 0;
        end
    end
    % Initialize nx3 matrix of all targeted non-zero neuron triplets.
    for i = length_vector
        if size(unique(input_timeseries(:,i)),1) == 1
            length_vector(i) = 0;
        end
    end
    length_vector(length_vector==0) = [];
    target_1 = nchoosek(length_vector,3);
    target_2 = circshift(target_1,1,2);
    target_3 = circshift(target_1,-1,2);
    all_triplets_list = [target_1; target_2; target_3];
    % Initialize output whose rows contain neuron triplets ijk that satisfy
    % TE(j->i)>TE(i->j) and TE(k->i)>TE(i->k) at the given time-delay.
    % Columns indicate in increasing order:
    % target | source1 | source2
    syms n
    functional_triplets = zeros(double(symsum(nchoosek(n,2), n, 2, size(length_vector,2)-1)), 3);
    % Initialize dummy variable to indicate rows of output matrix.
    row_index = 1;
    for i = 1:(size(all_triplets_list,1))
        if (threshold_matrix(all_triplets_list(i,1), all_triplets_list(i,2))>0) && (threshold_matrix(all_triplets_list(i,1), all_triplets_list(i,3))>0)
            functional_triplets(row_index,:) = all_triplets_list(i,:); % Write to row of output matrix indicated by row_index.
            row_index = row_index+1; % Increment dummy variable by 1.
        end
    end
    functional_triplets(functional_triplets(:,1)==0,:) = []; % Remove zeroes.
end