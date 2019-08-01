% Given two time-series T, S and one integer time-delay in units of time points, this function
% returns the transfer entropy from S to T with the given time-delay
% between T_future and T_past. Transfer entropy normalized by entropy of
% the target time-series is also calculated.
%
% Output is a scalar in units of bits.
%
% This function can only take discrete, scalar-valued time-series.

function [transfer_entropy, normed_TE] = TE(target, source, delay)
     % Check if inputs are of acceptable type and length.
    if (~isvector(target)) || (~isvector(source))
        error('Input time-series is not a vector.')
    elseif length(target) ~= length(source)
        error('Time-series are not of equal length.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay is not a positive integer.')
    end
    % Ensure time-series are represented as column vectors.
    target = target(:);
    source = source(:);
    % Truncate at beginning or end of time-series to create time-delay.
    target_future = target;
    target_future(1:delay) = [];
    target_past = target;
    target_past((length(target)-delay+1):length(target)) = [];
    source_past = source;
    source_past((length(source)-delay+1):length(source)) = [];
    transfer_entropy = cond_MI(source_past, target_future, target_past);
    % Normalize by entropy of target time-series.
    target_entropy = 0;
    for i = unique(target_future')
        prob = sum(target_future==i)/length(target_future);
        target_entropy = target_entropy - prob * log(prob) / log(2);
    end
    if target_entropy == 0
        disp('Target time-series has zero entropy. Using unnormalized transfer entropy.')
        normed_TE = transfer_entropy;
    else
        normed_TE = transfer_entropy / target_entropy;
    end
end