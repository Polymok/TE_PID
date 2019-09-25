% Timme et al., High-Degree Neurons Feed Cortical Computations (PLoS Computational Biology, 2016).
% 
% This function takes as input a target time-series T, a scalar value t
% that the target time-series may obtain, a time-delay between T_future and
% T_past, as well as an optional source time-series S, and returns the
% specific information, denoted I_spec(T_future=t; T_past, S).
%
% Output is a scalar in units of bits.
%
% This function can only take discrete, binary time-series.

function [spec_info] = I_spec(target, target_future_value, delay, opt_source)

    %% Check inputs.
    if ~isscalar(target_future_value)
        error('Input specific value is not a scalar.')
    elseif ~isscalar(delay)
        error('Input time-delay is not a scalar.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay is not a positive integer.')
    end 
    % Check if optional argument is given.
    if nargin == 4
        % Check if time-series inputs are of acceptable type and length.
        if (~isvector(target)) || (~isvector(opt_source))
            error('Input time-series are not vectors.')
        elseif (length(target)~=length(opt_source))
            error('Input time-series are not of the same length.')
        end
    elseif nargin == 3
        opt_source = zeros(length(target),1);
    else
        error('Number of arguments must be 3 or 4.')
    end
    % Ensure time-series are column vectors.
    target = target(:);
    opt_source = opt_source(:);
    
    %% Create past and future time-series using given time-delay.
    target_future = target;
    target_future(1:delay) = [];
    target_past = target;
    target_past((length(target)-delay+1):length(target)) = [];
    clear target;
    opt_source((length(opt_source)-delay+1):length(opt_source)) = [];
    
    %% Calculate specific information.
    spec_info = 0;
    % Sum over all possible values that the two time-series may take.
    for i = unique(target_past, 'rows')'
        for j = unique(opt_source, 'rows')'
            % Count number of instances instead of probability.
            jointprob = sum(sum([opt_source target_past target_future]==[j' i' target_future_value],2)==3);
            futureprob = sum(target_future==target_future_value);
            pastprob = sum(sum([opt_source target_past]==[j' i'],2)==2);
            % Discard cases with probability zero.
            if jointprob == 0
                disp(['Pr(source,target_future,target_past=', num2str(j'), ',', num2str(target_future_value), ',', num2str(i'), ') is zero. Case discarded.'])
            elseif futureprob == 0
                disp(['Pr(target_future=', num2str(target_future_value), ') is zero. Case discarded.'])
            elseif pastprob == 0
                disp(['Pr(target_past=', num2str(i'), ') is zero. Case discarded.'])
            else
                spec_info = spec_info + jointprob / futureprob * log2(jointprob / pastprob / futureprob * length(target_future));
            end
        end
    end
    
end