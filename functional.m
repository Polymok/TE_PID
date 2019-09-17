% Given a matrix whose columns contain each time-series and a time-delay,
% this function returns a 3-column matrix whose columns indicate neuron
% indexes in increasing order:
% target | source1 | source2 
% satisfying TE(source1->target)>TE(target->source1) and
% TE(source2->target)>TE(target->source2) at the given time-delay.
%
% The functional network is calculated using pairwise transfer entropy with
% the given delay. Bidirectional connections are only possible if transfer
% entropy in both directions are equal and non-zero. An array of percentage
% differences between non-zero functional weights when swapping target and
% source is also returned. This percentage is taken over the larger of the
% two transfer entropy values. The number of active functional neurons, the
% number of functional connections, and the number of
% bidirectional connections are printed to the command window.
%
% Optionally input a threshold percentage. E.g. if input threshold is 80%,
% result will discard bottom 80% of non-zero functional connections.
%
% This function is designed to filter neuron triplets for transfer entropy
% partial information decomposition calculation.

function [functional_triplets, functional_matrix, weight_diff] = functional(input_timeseries, delay, threshold)
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
        clear str
    end
    % Initialize weight matrix of zeroes whose rows are source neuron
    % indices and columns are target neuron indices. For a given neuron
    % pair ij and a time-delay, if the greater transfer entropy obtains
    % when i lags behind j, we record TE(i->j) in element (i,j) of the
    % directional matrix. Conversely, if greater transfer entropy obtains
    % when j lags behind i, we record TE(j->i in element (j,i) of the
    % directional matrix.
    functional_matrix = zeros(size(input_timeseries,2), size(input_timeseries,2));
    % Initialize list of all undirected active neuron pairs.
    length_vector = 1:size(input_timeseries,2);
    for i = length_vector
        if size(unique(input_timeseries(:,i)),1) == 1
            length_vector(length_vector==i) = 0;
        end
    end
    length_vector(length_vector==0) = [];
    neuron_pairs = nchoosek(length_vector,2);
    % Initialize vector of differences between weights w(ij) and w(ji).
    weight_diff = zeros(size(neuron_pairs,1),1);
    izeros = false(size(weight_diff,1),1);
    for i = 1:size(neuron_pairs,1)
        functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) = TE(input_timeseries(:,neuron_pairs(i,2)), input_timeseries(:,neuron_pairs(i,1)), delay);
        functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1)) = TE(input_timeseries(:,neuron_pairs(i,1)), input_timeseries(:,neuron_pairs(i,2)), delay);
        if (functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))==0) && (functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))==0)
            izeros(i) = 1;
        else
            if functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) > functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))
                weight_diff(i) = (functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))-functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))) / functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2));
                functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1)) = 0;
            elseif functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) < functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))
                weight_diff(i) = (functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))-functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))) / functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1));
                functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) = 0;
            end
        end
    end
    clear input_timeseries
    % Threshold.
    if nargin == 3
        if ~isscalar(threshold)
            error('Threshold must be a scalar.')
        elseif (threshold>1) || (threshold<0)
            error('Threshold must be between 0 and 1.')
        else
            ordered_weights = sort(functional_matrix(:));
            ordered_weights(ordered_weights==0) = [];
            threshold_value = ordered_weights(floor(threshold*size(ordered_weights,1)));
            ithreshold = false(size(weight_diff,1),1);
            for i = 1:size(neuron_pairs,1)
                if ((functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))>0) && (functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))<threshold_value)) || ((functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))>0) && (functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))<threshold_value))
                    ithreshold(i) = 1;
                end
            end
            weight_diff(ithreshold|izeros) = [];
            functional_matrix(functional_matrix<threshold_value) = 0;
        end
        clear threshold_value
        clear ordered_weights
        clear izeros
        clear ithreshold
    else
        weight_diff(izeros) = [];
        clear izeros
    end
    num_active = sum(sum(functional_matrix>0)>0);
    num_connect = sum(sum(functional_matrix>0));
    num_bidirect = sum(sum(functional_matrix==functional_matrix'&functional_matrix>0))/2;
    disp(['Number of active functional neurons: ', num2str(num_active)]);
    disp(['Number of functional connections: ', num2str(num_connect)]);
    disp(['Number of bidirectional functional connections: ', num2str(num_bidirect)]);
    % Initialize nx3 matrix of all targeted non-zero neuron triplets, where column one indicates the target neuron.
    target_1 = nchoosek(length_vector,3);
    target_2 = circshift(target_1,1,2);
    target_3 = circshift(target_1,-1,2);
    all_triplets_list = [target_1; target_2; target_3];
    clear target_1
    clear target_2
    clear target_3
    % Initialize output whose rows contain neuron triplets ijk that satisfy
    % TE(j->i)>TE(i->j) and TE(k->i)>TE(i->k) at the given time-delay.
    % Columns indicate in increasing order:
    % target | source1 | source2
    syms n
    functional_triplets = zeros(double(symsum(nchoosek(n,2), n, 2, size(length_vector,2)-1)), 3);
    clear length_vector
    % Initialize dummy variable to indicate rows of output matrix.
    row_index = 1;
    for i = all_triplets_list'
        if (functional_matrix(i(2), i(1))>0) && (functional_matrix(i(3), i(1))>0)
            functional_triplets(row_index,:) = i'; % Write to row of output matrix indicated by row_index.
            row_index = row_index+1;
        end
    end
    functional_triplets(functional_triplets(:,1)==0,:) = []; % Remove extra rows.
end