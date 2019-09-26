% Timme et al., "High-Degree Neurons Feed Cortical Computations". PLoS Comput Biol (2016).
%
% This function admits as input matrices, or cells of matrices. In the
% latter case, the cell must be one-dimensional. Each cell is interpreted
% as containing time-series of all neurons recorded in a single trial. In
% both cases, each time-series must be contained in a column of the input
% matrices, i.e. columns represent neurons and rows represent observations
% at certain times. Each spike train for each neuron must be scalar-valued.
%
% A positive integer time-delay is required to compute transfer entropy
% values.
%
% Optionally, input a list of neuron triplet indices to only calculate PID
% for a subset of total possible triplets. This list must take the form of
% an nx3 matrix, where the first column represents the target neuron index.
% If the input dataset is a cell of matrices, the neuron triplet list must
% also be a cell of identical dimensions, each containing a nx3 matrix. If
% this optional argument is not given, this function will calculate PID for
% all possible triplets.
%
% Optionally, enter a positive integer time-resolution at which to time bin
% input data. Note that the time-delay used to calculate transfer entropy
% will apply to the time-series after time binning.
%
% This function writes to output file a list of unique, synergistic, and
% redundant partial information terms of transfer entropy between all
% possible neuron triplets. The redundant partial information is taken to
% be equal to the minimum information function as defined by Williams and
% Beer, 2010, modified in Timme et al., 2016.
%
% The columns of the output file indicate in increasing order:
% target | source1 | source2 | synergy | redundancy | unique1 | unique2
%
% Entropy for each neuron is also calculated and recorded using the
% neuron's spike train with the first term(s) removed according to the
% time-delay given.
%
% PID is calculated for all neuron triplets. If N simultaneously recorded
% time-series are given, N*(N-1)*(N-2)/2 possible triplets are returned.

function TE_PID(output_filename, input_timeseries, delay, triplet_list, resolution)

    %% Check inputs.
    if ~isscalar(delay)
        error('Input time-delay must be a scalar.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay must be a positive integer.')
    elseif ~isstring(output_filename) && ~ischar(output_filename)
        error('Output filename must be entered as a string or char.')
    end
    
    %% Case: input is a matrix containing a single trial.
    if ismatrix(input_timeseries)
        % Check if a single time-series is contained in a column.
        if size(input_timeseries,1) < size(input_timeseries,2)
            str = input('Input matrix has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
            if str == 'y'
                input_timeseries = input_timeseries';
            end
            clear str
        end
        
        %% Optionally, time-bin at given resolution.
        if nargin==5
            input_timeseries = timebin(input_timeseries, resolution);
            clear time_resolution
        end
        
        %% Write outputs to separate file.
        output_file = fopen(output_filename, 'a');
        
        %% Calculate neuron entropies, and remove inactive neurons from subsequent PID calculations.
        nNeuron = size(input_timeseries,2);
        % If optional triplet list is given, extract list of unique neuron indices.
        if (nargin > 3) && ~isempty(triplet_list)
            if ~ismatrix(triplet_list)
                error('List of neuron triplets must be a matrix.')
            elseif size(triplet_list,2) ~= 3
                error('List of neuron triplets must have 3 columns.')
            elseif any(unique(triplet_list) > nNeuron)
                error('Neuron indices in given triplet list must not be greater than total number of neurons in input dataset.')
            else
                neuron_list = (unique(triplet_list))';
            end
        % If optional triplet list is not given, create list of all neuron indices.
        else
            neuron_list = 1:nNeuron;
        end
        % Calculate and record entropy of target neuron separately.
        fprintf(output_file, 'Target, Entropy\n');
        for i = neuron_list
            entropy = 0;
            target_future = input_timeseries(:,i);
            target_future(1:delay) = [];
            for j = unique(target_future)'
                prob = sum(target_future==j')/size(target_future,1);
                if prob ~= 0
                    entropy = entropy - prob * log2(prob);
                end
            end
            fprintf(output_file, '%1u, %.6g\n', i, entropy);
            % Remove inactive neurons.
            if entropy == 0
                neuron_list(i) = 0;
            end
        end
        neuron_list(neuron_list==0) = [];
        clear prob entropy
        
        %% If optional triplet list is not given, create list of all active neuron triplets.
        if nargin==3 || isempty(triplet_list)
            target_1 = nchoosek(neuron_list,3);
            target_2 = circshift(target_1,1,2);
            target_3 = circshift(target_1,-1,2);
            triplet_list = [target_1; target_2; target_3];
            clear target_1 target_2 target_3;
        end
        
        %% Calculate and record PID values.
        fprintf(output_file, 'Target, Source1, Source2, Synergy, Redundancy, Unique1, Unique2\n');
        % Calculate and store transfer entropies from single source to single target.
        targeted_pairs = unique([triplet_list(:,1:2); triplet_list(:,1) triplet_list(:,3)], 'rows');
        single_TEs = zeros(size(targeted_pairs,1),3);
        row_index = 1;
        for i = targeted_pairs'
            single_TEs(row_index,1:2) = i';
            single_TEs(row_index,3) = TE(input_timeseries(:,i(1)), input_timeseries(:,i(2)), delay);
            row_index = row_index + 1;
        end
        clear neuron_list targeted_pairs row_index
        % Import all neuron triplets from triplet_list, and calculate PID values.
        for i = triplet_list'
            redundancy = I_min_TE(input_timeseries(:,i(1)), input_timeseries(:,i(2)), input_timeseries(:,i(3)), delay);
            TE1 = single_TEs(sum(single_TEs(:,1:2)==[i(1) i(2)],2)==2,3);
            TE2 = single_TEs(sum(single_TEs(:,1:2)==[i(1) i(3)],2)==2,3);
            TE12 = TE(input_timeseries(:,i(1)), [input_timeseries(:,i(2)) input_timeseries(:,i(3))], delay);
            unique1 = TE1 - redundancy;
            unique2 = TE2 - redundancy;
            synergy =  TE12 - redundancy - unique1 - unique2;
            output_data = [i; synergy; redundancy; unique1; unique2];
            fprintf(output_file, '%1u, %1u, %1u, %.6g, %.6g, %.6g, %.6g\n', output_data);
        end
        fclose(output_file);
        
    %% Case: input is a cell containing multiple trials.
    % Extract each cell element as a matrix, and feed back into TE_PID.m.
    elseif iscell(input_timeseries)
        if (size(input_timeseries,1) ~= 1) && (size(input_timeseries,2) ~= 1)
            error('Input cell must be one-dimensional.')
        elseif size(input_timeseries,1) > size(input_timeseries,2)
            input_timeseries = input_timeseries';
        end
        for trial = 1:size(input_timeseries,2)
            input_matrix = input_timeseries{1, trial};
            % Check if optional triplet list is given.
            if nargin == 3 || isempty(triplet_list)
                if nargin==5
                    TE_PID(output_filename, input_matrix, delay, [], resolution);
                else
                    TE_PID(output_filename, input_matrix, delay);
                end
            else
                % Check if given triplet list is a matrix or a cell.
                if ismatrix(triplet_list)
                    if nargin==5
                        TE_PID(output_filename, input_matrix, delay, triplet_list, resolution);
                    else
                        TE_PID(output_filename, input_matrix, delay, triplet_list);
                    end
                elseif iscell(triplet_list)
                    if (size(triplet_list,1) ~= 1) && (size(triplet_list,2) ~= 1)
                        error('Input cell of triplet lists must be one-dimensional.')
                    elseif size(triplet_list,1) > size(triplet_list,2)
                        triplet_list = triplet_list';
                    end
                    triplet_matrix = triplet_list{1, trial};
                    if nargin==5
                        TE_PID(output_filename, input_matrix, delay, triplet_matrix, resolution);
                    else
                        TE_PID(output_filename, input_matrix, delay, triplet_matrix);
                    end
                end
            end
        end
    else
        error('Input dataset must be a cell or a matrix.')
    end
    
end