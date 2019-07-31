% Given two time-series and a maximum time delay constraint, this function
% returns the time lag between the two time-series smaller than the maximal
% constraint that yields the greatest transfer entropy normalized by the
% target time-series's entropy.
%
% This function can only take scalar-valued time-series.

function [timelagged_TE_normed, timelag] = TE_timelag(u,v,range)
    if (~isvector(u)) || (~isvector(v))
        error('Input is not a vector.')
    elseif ~isscalar(range)
        error('Maximum time lag is not a scalar.')
    elseif length(u) ~= length(v)
        error('Inputs not of the same length.')
    end
    u = u(:); % Ensure inputs are column vectors.
    v = v(:);
    normed_transfers = zeros(1, 2*floor(range));
    for i = 1:floor(range)
        [~, normed_transfers(i)] = TE(u,v,i);
        [~, normed_transfers(i+floor(range))] = TE(v,u,i);
    end
    [timelagged_TE_normed, timelag] = max(normed_transfers);
    if timelag > floor(range)
        timelag = floor(range)-timelag;
    end
end
