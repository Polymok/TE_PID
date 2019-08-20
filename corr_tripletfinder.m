% Given a matrix of time-series and a time-delay, this function returns an
% asymmetric weight matrix of correlation coefficients between each pair of
% time-series. For each time-series pair ij, if the greater correlation
% coefficient obtains when j lags behind i, matrix element (i,j) is set to
% 0. If i lags behind j, matrix element (j,i) is set to 0.
%
% This function further returns a list of time-series triplets ijk that
% satisfies weight_matrix(j,i)>0 and weight_matrix(k,i)>0.
%
% Pearson's correlation coefficient is used to construct the weight matrix.
%
% This function can only take scalar-valued time-series.

function [output_triplets, weight_matrix] = corr_tripletfinder(input_timeseries, delay)
    if ~ismatrix(input_timeseries)
        error('Input is not a matrix.')
    elseif ~isscalar(delay)
        error('Time-delay is not a scalar.')
    elseif size(input_timeseries,1) < size(input_timeseries,2)
        str = input('Input matrix has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            input_timeseries = input_timeseries';
        end
    end
    weight_matrix = corrcoef(input_timeseries);
    for i = 1:size(input_timeseries,2)
        for j = (i+1):size(input_timeseries,2)
            i_past = input_timeseries(:,i);
            i_past((size(input_timeseries,1)-delay+1):size(input_timeseries,1),:) = [];
            i_future = input_timeseries(:,i);
            i_future(1:delay,:) = [];
            j_past = input_timeseries(:,j);
            j_past((size(input_timeseries,1)-delay+1):size(input_timeseries,1),:) = [];
            j_future = input_timeseries(:,j);
            j_future(1:delay,:) = [];
            i_to_j = corrcoef(i_past, j_future);
            j_to_i = corrcoef(j_past, i_future);
            if i_to_j(1,2) > j_to_i(1,2)
                weight_matrix(j,i) = 0;
            elseif i_to_j(1,2) < j_to_i(1,2)
                weight_matrix(i,j) = 0;
            end
        end
    end
    output_triplets = zeros(size(input_timeseries,2)*(size(input_timeseries,2)-1)*(size(input_timeseries,2)-2)/4, 3);
    row_index = 1;
    for i = 1:(size(input_timeseries,2)) % Target neuron.
        for j = 1:(size(input_timeseries,2)-1) % Source neuron 1.
            if i==j
            else
                for k = (j+1):(size(input_timeseries,2)) % Source neuron 2.
                    if i==k
                    else
                        if (weight_matrix(j,i)>0) && (weight_matrix(k,i)>0)
                            output_triplets(row_index,:) = [i j k];
                            row_index = row_index+1;
                        end
                    end
                end
            end
        end
    end
end
