% Given a matrix whose columns contain each time-series and a maximum
% time-delay, this function returns a 3-column matrix whose columns
% indicate neuron indexes in increasing order: 
% target | source1 | source2
% satisfying TE(source1->target)>TE(target->source1) and
% TE(source2->target)>TE(target->source2) at all time-delays less than the
% maximal time-delay constraint.
%
% This function is designed to filter neuron triplets for transfer entropy
% partial information decomposition calculation.

function [output_triplets] = TE_tripletfinder(input_timeseries, max_delay)
    % Check if input formats are acceptable.
    if ~isa(input_timeseries, 'double') && ~isa(input_timeseries, 'single') && ~isa(input_timeseries, 'logical')
        error('Input time-series must be a matrix.')
    elseif ~isscalar(max_delay)
        error('Input maximum time-delay must be a scalar.')
    elseif max_delay<1
        error('Input maximum time-delay must be at least 1.')
    end
    % Check if a single time-series is contained in a column.
    if size(input_timeseries,1) < size(input_timeseries,2)
        str = input('Input matrix has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            input_timeseries = input_timeseries';
        end
        clear str
    end
    % Initialize directional matrix of zeroes whose rows are source neuron
    % indices and columns are target neuron indices. For a given neuron
    % pair ij and a maximal time-delay, if maximal transfer entropy obtains
    % when i lags behind j, we record a 1 in element (i,j) of the
    % directional matrix. Conversely, if maximal transfer entropy obtains
    % when j lags behind i, we record a 1 in element (j,i) of the
    % directional matrix.
    directional_matrix = zeros(size(input_timeseries,2), size(input_timeseries,2));
    for i = 1:(size(input_timeseries,2)-1)
        for j = (i+1):(size(input_timeseries,2))
            [~, direction_ij] = TE_timelag(input_timeseries(:,i), input_timeseries(:,j), max_delay);
            if direction_ij > 0
                directional_matrix(j,i) = 1;
            elseif direction_ij < 0
                directional_matrix(i,j) = 1;
            end
        end
    end
    % Initialize output whose rows contain neuron triplets ijk that satisfy
    % TE(j->i)>TE(i->j) and TE(k->i)>TE(i->k) at time-delays less than the
    % maximal constraint. Columns indicate in increasing order:
    % target | source1 | source2
    output_triplets = zeros(size(input_timeseries,2)*(size(input_timeseries,2)-1)*(size(input_timeseries,2)-2)/4, 3);
    row_index = 1;
    for i = 1:(size(input_timeseries,2)) % Target neuron.
        for j = 1:(size(input_timeseries,2)-1) % Source neuron 1.
            if i==j
            else
                for k = (j+1):(size(input_timeseries,2)) % Source neuron 2.
                    if i==k
                    else
                        if (directional_matrix(j,i)==1) && (directional_matrix(k,i)==1)
                            output_triplets(row_index,:) = [i j k];
                            row_index = row_index+1;
                        end
                    end
                end
            end
        end
    end    
end