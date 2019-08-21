% Timme et al., "High-Degree Neurons Feed Cortical Computations". PLoS Comput Biol (2016).
%
% This function admits as inputs matrices, or cells of matrices. In the
% latter case, the cell must be one-dimensional. Each cell is interpreted
% as containing time-series of all neurons recorded in a single trial. In
% both cases, each time-series must be contained in a column of the input
% matrices, i.e. columns represent neurons and rows represent observations
% at certain times.
%
% A positive integer time-delay is required to compute transfer entropy
% values.
%
% Optionally, you may input a list of neuron triplet indices if you wish to
% only calculate PID for a subset of total possible triplets. This list
% must take the form of an nx3 matrix, where the first column represents
% the target neuron index. If the input dataset is a cell of matrices, the
% neuron triplet list must also be a cell of identical dimensions, each
% containing a nx3 matrix. If this optional argument is not given, this
% function will calculate PID for all possible triplets.
%
% This function returns as output a matrix of unique, synergistic, and
% redundant partial information terms of transfer entropy between all
% possible neuron triads. The redundant partial information is taken to be
% equal to the minimum information function as defined by Williams and
% Beer, 2010, modified in Timme et al., 2016.
%
% In each output matrix, the columns indicate in increasing order:
% target | source1 | source2 | synergy | redundancy | unique1 | unique2
%
% PID is calculated for all neuron triplets. If N simultaneously recorded
% time-series are given, N*(N-1)*(N-2)/2 triplets are returned.

function [output_data] = TE_PID(input_data, delay, triplet_list)
    if ~isscalar(delay)
        error('Input time-delay must be a scalar.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay must be a positive integer.')
    end
    % Case: input is a matrix containing a single trial.
    if ismatrix(input_data)
        % Check if a single time-series is contained in a column.
        if size(input_data,1) < size(input_data,2)
            str = input('Input matrix has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
            if str == 'y'
                input_data = input_data';
            end
        end
        % Check if optional triplet list is given.
        if nargin == 2 || isempty(triplet_list)
            vector_neurons = 1:size(input_data,2);
            target_1 = nchoosek(vector_neurons,3);
            target_2 = circshift(target_1,1,2);
            target_3 = circshift(target_1,-1,2);
            triplet_list = [target_1; target_2; target_3];
        elseif nargin == 3
            if ~ismatrix(triplet_list)
                error('List of neuron triplets must be a matrix.')
            elseif size(triplet_list,2) ~= 3
                error('List of neuron triplets must have 3 columns.')
            elseif any(unique(triplet_list) > size(input_data,2))
                error('Neuron indices in given triplet list must not be greater than total number of neurons in input dataset.')
            end
        end
        % Time bin functionality.
%         str = input('Time bin input time-series? y/n: ', 's');
%         if str == 'y'
%             resolution = input('Choose a time resolution. Please input a positive integer: ');
%             while ~isscalar(resolution)
%                 resolution = input('Time resolution must be a scalar. Please input a positive integer: ');
%             end
%             input_data = timebin(input_data, resolution);
%         end
%         clear str
        output_data = zeros(size(triplet_list,1), 7); % Initialize output matrix for a single trial. 7 total indices: target, source1, source2, synergy, redundancy, unique1, unique2.
        % Initialize row index to write outputs to output_matrix.
        row_index = 1;
        % Import all neuron triplets from triplet_list.
        for i = 1:(size(triplet_list,1))
            redundancy = I_min_TE(input_data(:,triplet_list(i,1)), input_data(:,triplet_list(i,2)), input_data(:,triplet_list(i,2)), delay);
            [~, TE1] = TE(input_data(:,triplet_list(i,1)), input_data(:,triplet_list(i,2)), delay); % Use transfer entropy normalized by entropy of target.
            [~, TE2] = TE(input_data(:,triplet_list(i,1)), input_data(:,triplet_list(i,3)), delay);
            [~, TE12] = TE(input_data(:,triplet_list(i,1)), [input_data(:,triplet_list(i,2)) input_data(:,triplet_list(i,3))], delay);
            unique1 = TE1 - redundancy;
            unique2 = TE2 - redundancy;
            synergy =  TE12 - redundancy - unique1 - unique2;
            output_data(row_index,:) = [triplet_list(i,:) synergy redundancy unique1 unique2]; % Write to row of output_matrix indicated by row_index.
            row_index = row_index+1; % Increment row index by 1.
        end
    % Case: input is a cell containing multiple trials.
    elseif iscell(input_data)
        if (size(input_data,1) ~= 1) && (size(input_data,2) ~= 1)
            error('Input cell must be one-dimensional.')
        elseif size(input_data,1) > size(input_data,2)
            input_data = input_data';
        end
        output_data = cell(1, size(input_data,2)); % Initialize output cell array.
        for trial = 1:size(input_data,2)
            input_matrix = input_data{1, trial};
            % Check if optional triplet list is given.
            if nargin == 2
                output_matrix = TE_PID(input_matrix, delay); % Feed extracted matrix back into function.
            elseif nargin == 3
                % Check if given triplet list is a matrix or a cell.
                if ismatrix(triplet_list)
                    output_matrix = TE_PID(input_matrix, delay, triplet_list);
                elseif iscell(triplet_list)
                    if (size(triplet_list,1) ~= 1) && (size(triplet_list,2) ~= 1)
                        error('Input cell of triplet lists must be one-dimensional.')
                    elseif size(triplet_list,1) > size(triplet_list,2)
                        triplet_list = triplet_list';
                    end
                    triplet_matrix = triplet_list{1, trial};
                    output_matrix = TE_PID(input_matrix, delay, triplet_matrix);
                end
            end
            output_data{1, trial} = output_matrix; % Write resulting PID values for each trial back into output cell.
        end
    else
        error('Input dataset must be a cell or a matrix.')
    end
end