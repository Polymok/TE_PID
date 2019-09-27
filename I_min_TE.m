% Timme et al., High-Degree Neurons Feed Cortical Computations (PLoS Computational Biology, 2016).
% 
% This function takes as input two source time-series S_1 and S_2, one
% target time-series T, and a time-delay.
%
% This function returns the minimum information from the two sources to the
% target future conditioned on the target past, denoted 
% I_min(T_F; S_1, S_2 | T_P).
%
% Output is a scalar in units of bits.
%
% This function can only take discrete time-series.

function [min_info] = I_min_TE(target, source1, source2, delay)

    %% Check inputs.
    if nargin < 4
        error('4 inputs required.')
    elseif (~isvector(target)) || (~isvector(source1)) || (~isvector(source2))
        error('Inputs are not vectors.')
    elseif (length(target)~=length(source1)) || (length(target)~=length(source2))
        error('Input vectors are not of equal length.')
    elseif ~isscalar(delay)
        error('Input time-delay is not a scalar.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay is not a positive integer.')
    end
    % Ensure inputs are column vectors.
    target = target(:);
    source1 = source1(:);
    source2 = source2(:);
    
    %% Create future target time-series using given time-delay.
    target_future = target;
    target_future(1:delay) = [];
    
    %% Calculate minimum information.
    min_info = 0;
    for i = unique(target_future, 'rows')'
        subtractor = I_spec(target, i', delay);
        spec_info1 = I_spec(target, i', delay, source1) - subtractor;
        spec_info2 = I_spec(target, i', delay, source2) - subtractor;
        min_info = min_info + sum(target_future==i')/length(target_future)*min(spec_info1, spec_info2);
    end
    
end