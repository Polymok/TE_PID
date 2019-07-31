% Given a 1-dimensional time-series T, a 2-dimensional time-series S, and
% one integer time-delay, this function returns the transfer entropy from S
% to T with the given time-delay between T_future and T_past. Transfer
% entropy normalized by entropy of the target time-series is also
% calculated.
%
% Output is a scalar in units of bits.
%
% This function can only take discrete time-series. The first input
% time-series must be 1-dimensional. The second input must be
% 2-dimensional. This design is particular to calculating the transfer
% entropy partial information decomposition.

function [transfer_entropy, normed_TE] = TE_2dim(target, source, delay)
     % Check if inputs are of acceptable type and length.
    if (~isvector(target))
        error('Input target time-series is not 1-dimensional.')
    elseif (size(source,1)~=2) && (size(source,2)~=2)
        error('Input source time-series is not 2-dimensional.')
    elseif length(target) ~= length(source)
        error('Time-series are not of equal length.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay is not a positive integer.')
    end
    % Ensure time-series are represented as column vectors.
    target = target(:);
    if size(source,1) < size(source,2)
        str = input('Input 2-dimensional time-series has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input time-series? y/n: ','s');
        if str == 'y'
            source = source';
        end
        clear str
    end
    % Remove elements at beginning or end of time-series to create time-delay.
    target_future = target;
    target_future(1:delay) = [];
    target_past = target;
    target_past((length(target)-delay+1):length(target)) = [];
    source1_past = source(:,1);
    source1_past((length(source1_past)-delay+1):length(source1_past)) = [];
    source2_past = source(:,2);
    source2_past((length(source2_past)-delay+1):length(source2_past)) = [];
    transfer_entropy = cond_MI_2dim([source1_past source2_past], target_future, target_past);
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