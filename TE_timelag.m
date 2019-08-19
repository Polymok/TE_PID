% Given two time-series u,v and a maximum time-delay constraint, this
% function returns the time-lag between the two time-series smaller than
% the maximum time-delay constraint that yields the greatest transfer
% entropy normalized by the target time-series's entropy.
% 
% A positive time-lag value means maximal transfer entropy obtains when v
% lags behind u. A negative time-lag value means maximal transfer entropy
% obtains when u lags behind v. The structure of transfer entropy procludes
% time-lag values of 0.
%
% This function can only take discrete time-series.

function [timelagged_TE_normed, timelag] = TE_timelag(u,v,max_delay)
    % Ensure inputs are column vectors.
    if size(u,1) < size(u,2)
        str = input('Input vector has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            u = u';
        end
    end
    if size(v,1) < size(v,2)
        str = input('Input vector has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            v = v';
        end
    end
    if ~isscalar(max_delay)
        error('Maximum time-lag is not a scalar.')
    elseif max_delay<1
        error('Maximum time-lag must be at least 1.')
    elseif size(u,1) ~= size(v,1)
        error('Inputs not of the same length.')
    end
    normed_transfers = zeros(1, 2*floor(max_delay));
    for i = 1:floor(max_delay)
        [~, normed_transfers(i)] = TE(u,v,i);
        [~, normed_transfers(i+floor(max_delay))] = TE(v,u,i);
    end
    [timelagged_TE_normed, timelag] = max(normed_transfers);
    if timelag > floor(max_delay)
        timelag = floor(max_delay)-timelag;
    end
end
