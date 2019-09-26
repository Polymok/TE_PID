% Given a matrix whose columns contain each time-series and a time-delay,
% return a 3-column matrix whose columns indicate neuron indexes:
% target | source1 | source2 
% satisfying TE(source1->target)>TE(target->source1) and
% TE(source2->target)>TE(target->source2) at the given time-delay.
%
% Calculate functional network using pairwise transfer entropy with the
% given delay. Bidirectional connections are only possible if transfer
% entropy in both directions are equal and non-zero. A list of percentage
% differences between non-zero functional weights when swapping target and
% source is also returned. This percentage is taken over the smaller of the
% two transfer entropy values. The number of active functional neurons, the
% number of functional connections, and the number of bidirectional
% connections are printed to the command window.
%
% Optionally input a threshold percentage. E.g. if input threshold is 80%,
% result will discard bottom 80% of non-zero functional connections.
%
% This function is designed to filter neuron triplets for transfer entropy
% partial information decomposition calculation.

function [functional_triplets, functional_matrix, weight_diff, raw_TE_matrix] = functional(input_timeseries, delay, threshold)

    %% Check inputs.
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
    nNeuron = size(input_timeseries,2);
    
    %% Calculate functional weight matrix using transfer entropy.
    % Initialize weight matrix whose rows are source neuron
    % indices and columns are target neuron indices. For a given neuron
    % pair ij and a time-delay, if the greater transfer entropy obtains
    % when i lags behind j, record TE(i->j) in element (i,j) of the
    % directional matrix. Conversely, if greater transfer entropy obtains
    % when j lags behind i, we record TE(j->i in element (j,i) of the
    % directional matrix.
    functional_matrix = zeros(nNeuron, nNeuron);
    % Initialize list of all undirected active neuron pairs.
    neuron_list = 1:nNeuron;
    for i = neuron_list
        if size(unique(input_timeseries(:,i)),1) == 1
            neuron_list(neuron_list==i) = 0;
        end
    end
    neuron_list(neuron_list==0) = [];
    neuron_pairs = nchoosek(neuron_list,2);
    for i = 1:size(neuron_pairs,1)
        functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) = TE(input_timeseries(:,neuron_pairs(i,2)), input_timeseries(:,neuron_pairs(i,1)), delay);
        functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1)) = TE(input_timeseries(:,neuron_pairs(i,1)), input_timeseries(:,neuron_pairs(i,2)), delay);
    end
    clear input_timeseries
    
    %% Record non-directed weight matrix.
    if nargout == 4
        raw_TE_matrix = functional_matrix;
    end
    
    %% Record difference in weights when swapping target and source. Set smaller of the two to zero.
    weight_diff = zeros(size(neuron_pairs,1),1);
    % If either w(ij) or w(ji) is zero, exclude pair from weight difference list.
    iZeros = false(size(weight_diff,1),1);
    for i = 1:size(neuron_pairs,1)
        if (functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))==0) || (functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))==0)
            iZeros(i) = 1;
        elseif functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) > functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))
            weight_diff(i) = (functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))-functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))) / functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1));
            functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1)) = 0;
        elseif functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) < functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))
            weight_diff(i) = (functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))-functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))) / functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2));
            functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) = 0;
        elseif functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2)) == functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))
            weight_diff(i) = eps; % If weight pairs are non-zero and equal, record eps to facilitate taking log plots.
        end
    end
    
    %% Optionally, threshold functional weight matrix.
    if nargin == 3
        if ~isscalar(threshold)
            error('Threshold must be a scalar.')
        elseif (threshold>1) || (threshold<0)
            error('Threshold must be between 0 and 1.')
        else
            ordered_weights = sort(functional_matrix(:));
            ordered_weights(ordered_weights==0) = [];
            threshold_value = ordered_weights(floor(threshold*size(ordered_weights,1)));
            % Remove weight differences corresponding to thresholded functional connections.
            iThreshold = false(size(weight_diff,1),1);
            for i = 1:size(neuron_pairs,1)
                if ((functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))>0) && (functional_matrix(neuron_pairs(i,1),neuron_pairs(i,2))<threshold_value)) || ((functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))>0) && (functional_matrix(neuron_pairs(i,2),neuron_pairs(i,1))<threshold_value))
                    iThreshold(i) = 1;
                end
            end
            clear neuron_pairs
            weight_diff(iThreshold|iZeros) = [];
            functional_matrix(functional_matrix<threshold_value) = 0;
        end
        clear threshold_value ordered_weights izeros ithreshold
    else
        weight_diff(iZeros) = [];
        clear izeros
    end
    
    %% Print topological data of interest to command window.
    disp(['Number of active functional neurons: ', num2str(sum(sum(functional_matrix>0)>0))]);
    disp(['Number of functional connections: ', num2str(sum(sum(functional_matrix>0)))]);
    disp(['Number of bidirectional functional connections: ', num2str(sum(sum(functional_matrix==functional_matrix'&functional_matrix>0))/2)]);
    
    %% Compute functional triplets from weight matrix.
    % Initialize nx3 matrix of all targeted non-zero neuron triplets, where column one indicates the target neuron.
    target_1 = nchoosek(neuron_list,3);
    target_2 = circshift(target_1,1,2);
    target_3 = circshift(target_1,-1,2);
    all_triplets_list = [target_1; target_2; target_3];
    clear target_1 target_2 target_3
    % Initialize output whose rows contain neuron triplets ijk that satisfy
    % TE(j->i)>TE(i->j) and TE(k->i)>TE(i->k) at the given time-delay.
    % Columns indicate in increasing order:
    % target | source1 | source2
    functional_triplets = zeros(nNeuron*(nNeuron-1)*(nNeuron-2)/2, 3); % Upper bound on number of possible fan-in triplets.
    clear neuron_list
    % Initialize dummy variable to indicate rows of output matrix.
    row_index = 1;
    for i = all_triplets_list'
        if (functional_matrix(i(2),i(1))>0) && (functional_matrix(i(3),i(1))>0)
            functional_triplets(row_index,:) = i'; % Write to row of output matrix indicated by row_index.
            row_index = row_index+1;
        end
    end
    functional_triplets(functional_triplets(:,1)==0,:) = []; % Remove extra rows.
    
end