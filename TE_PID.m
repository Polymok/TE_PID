% Timme et al., "High-Degree Neurons Feed Cortical Computations". PLoS Comput Biol (2016).
%
% This function admits as inputs matrices, or cells of matrices. In the
% latter case, the cell must be one-dimensional. Each cell is interpreted
% as containing time-series of all neurons recorded in a single trial. In
% both cases, each time-series must be contained in a column of the input
% matrices.
%
% A positive integer time-delay is required to compute transfer entropy
% values.
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
%
% TODO:
% Add first round filtering of neuron triplets by identifying significant
% transfer entropy values.

function [output_data] = TE_PID(input_data, delay)
    if ~isscalar(delay)
        error('Input time-delay is not a scalar.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay is not a positive integer.')
    end
    % Case: input is a matrix containing a single trial.
    if isa(input_data, 'double') || isa(input_data, 'single')
        % Check if a single time-series is contained in a column.
        if size(input_data,1) < size(input_data,2)
            str = input_data('Input matrix has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
            if str == 'y'
                input_data = input_data';
            end
            clear str
        end
        output_data = zeros(size(input_data,2)*(size(input_data,2)-1)*(size(input_data,2)-2)/2, 7); % Initialize output matrix for a single trial. 7 total indices: target, source1, source2, synergy, redundancy, unique1, unique2.
        % Pick out all neuron triplets.
        row_index = 1; % Initialize row index to write outputs to output_matrix. 
        for i = 1:(size(input_data,2)) % Target neuron.
            for j = 1:(size(input_data,2)-1) % Source neuron 1.
                if i==j
                else
                    for k = (j+1):(size(input_data,2)) % Source neuron 2.
                        if i==k
                        else
                            redundancy = I_min_TE(input_data(:,i), input_data(:,j), input_data(:,k), delay);
                            [~, TE1] = TE(input_data(:,i), input_data(:,j), delay); % Use transfer entropy normalized by entropy of target.
                            [~, TE2] = TE(input_data(:,i), input_data(:,k), delay);
                            [~, TE12] = TE_2dim(input_data(:,i), [input_data(:,j) input_data(:,k)], delay);
                            unique1 = TE1 - redundancy;
                            unique2 = TE2 - redundancy;
                            synergy =  TE12 - redundancy - unique1 - unique2;
                            output_data(row_index,:) = [i j k synergy redundancy unique1 unique2]; % Write to row of output_matrix indicated by row index.
                            row_index = row_index+1; % Increment row index by 1.
                        end
                    end
                end
            end
        end
    % Case: input is a cell containing multiple trials.
    elseif isa(input_data, 'cell')
        if (size(input_data,1) ~= 1) && (size(input_data,2) ~= 1)
            error('Input is not a one-dimensional cell.')
        elseif size(input_data,1) > size(input_data,2)
        input_data = input_data';
        end
        output_data = cell(1, size(input_data,2)); % Initialize output cell array.
        for trial = 1:size(input_data,2)
            input_matrix = input_data{1, trial};
            % Check if a single time-series is contained in a column.
            if size(input_matrix,1) < size(input_matrix,2)
                str = input_data('Input matrix has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
                if str == 'y'
                    input_matrix = input_matrix';
                end
                clear str
            end
            output_matrix = zeros(size(input_matrix,2)*(size(input_matrix,2)-1)*(size(input_matrix,2)-2)/2, 7); % Initialize output matrix for a single trial. 7 total indices: target, source1, source2, synergy, redundancy, unique1, unique2.
            % Pick out all neuron triplets.
            row_index = 1; % Initialize row index to write outputs to output_matrix. 
            for i = 1:(size(input_matrix,2)) % Target neuron.
                for j = 1:(size(input_matrix,2)-1) % Source neuron 1.
                    if i==j
                    else
                        for k = (j+1):(size(input_matrix,2)) % Source neuron 2.
                            if i==k
                            else
                                redundancy = I_min_TE(input_matrix(:,i), input_matrix(:,j), input_matrix(:,k), delay);
                                [~, TE1] = TE(input_matrix(:,i), input_matrix(:,j), delay); % Use transfer entropy normalized by entropy of target.
                                [~, TE2] = TE(input_matrix(:,i), input_matrix(:,k), delay);
                                [~, TE12] = TE_2dim(input_matrix(:,i), [input_matrix(:,j) input_matrix(:,k)], delay);
                                unique1 = TE1 - redundancy;
                                unique2 = TE2 - redundancy;
                                synergy =  TE12 - redundancy - unique1 - unique2;
                                output_matrix(row_index,:) = [i j k synergy redundancy unique1 unique2]; % Write to row of output_matrix indicated by row index.
                                row_index = row_index+1; % Increment row index by 1.
                            end
                        end
                    end
                end
            end
            output_data{1, trial} = output_matrix;
        end
    else
    error('Input dataset is not a cell or matrix.')
    end
end