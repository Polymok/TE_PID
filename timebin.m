% Given a positive integer time resolution and a binary matrix whose
% columns contain the time-series for a single neuron, this function
% returns a binary matrix binned at the given temporal resolution.
%
% The input time-series are partitioned into lengths equal to the given
% time resolution, with the first partition containing the 1-st to the
% (time resolution+1)-th value of the input time-series. If the n-th
% partition contains at least a single 1 value, a 1 is written to the n-th
% row of the output time-series. Otherwise, a 0 is written.
%
% NOTE: This function truncates the end of the time-series if the length of
% the time-series is not a multiple of the time resolution. Namely, if a
% time-series of length N and a time resolution of T is given, the last (N
% mod T) input values are discarded.

function [output_timeseries] = timebin(input_timeseries, resolution)

%     %% Check inputs.
%     if ~isscalar(resolution)
%         error('Input time-delay is not a scalar.')
%     elseif (round(resolution)~=resolution) || (resolution<1)
%         error('Input time-delay is not a positive integer.')
%     end
%     % Check if input time-series is a matrix.
%     if ~ismatrix(input_timeseries)
%         error('Input time-series must be a matrix.')
%     end
%     % Check if time-series are binary.
%     entries = unique(input_timeseries);
%     if size(entries,1) ~= 2
%         error('Input time-series must be binary.')
%     elseif entries ~= [0;1]
%         str = input('Input matrix must be binary with 0 or 1 valued entries. Change matrix to contain only 0 and 1? y/n: ','s');
%         if str == 'y'
%             input_timeseries = input_timeseries==entries(2);
%         end
%     end
%     % Check if a single time-series is contained in a column.
%     if size(input_timeseries,1) < size(input_timeseries,2)
%         str = input('Input matrix has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
%         if str == 'y'
%             input_timeseries = input_timeseries';
%         end
%         clear str
%     end
    
    %% Time bin.
    % Initialize output to contain rows equal to the number of input rows divided by the resolution rounded down.
    output_timeseries = zeros(floor(size(input_timeseries,1)/resolution), size(input_timeseries,2));
    for i = 1:size(output_timeseries,1)
        for j = 1:size(input_timeseries,2)
            if sum(input_timeseries(((i-1)*resolution+1):(i*resolution),j)==1) ~= 0
                output_timeseries(i,j) = 1;
            else
                output_timeseries(i,j) = 0;
            end
        end
    end
    
end