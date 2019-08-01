% Timme et al., High-Degree Neurons Feed Cortical Computations (PLoS Computational Biology, 2016).
% 
% This function takes as input a target time-series T, a scalar value t
% that the target time-series may obtain, a time-delay between T_future and
% T_past, as well as an optional source time-series S, and returns the
% specific information, denoted I_spec(T_future=t; T_past, S).
%
% Output is a scalar in units of bits.
%
% This function can only take discrete, scalar-valued time-series.

function [spec_info] = I_spec(target, target_future_value, delay, opt_source)
    % Check if specific value and time-delay are scalars.
    if ~isscalar(target_future_value)
        error('Input specific value is not a scalar.')
    elseif ~isscalar(delay)
        error('Input time-delay is not a scalar.')
    elseif (round(delay)~=delay) || (delay<1)
        error('Input time-delay is not a positive integer.')
    end
    % Ensure target time-series is a column vector.
    target = target(:);
    % Initialize output.
    spec_info = 0; 
    % Check if optional argument is given.
    if nargin == 4
        % Check if time-series inputs are of acceptable type and length.
        if (~isvector(target)) || (~isvector(opt_source))
            error('Input time-series are not vectors.')
        elseif (length(target)~=length(opt_source))
            error('Input time-series are not of the same length.')
        end
        % Ensure optional time-series is a column vector.
        opt_source = opt_source(:);
        % Create past and future time-series for target using given
        % time-delay by removing elements at start or end of time-series.
        target_future = target;
        target_future(1:delay) = [];
        target_past = target;
        target_past((length(target)-delay+1):length(target)) = [];
        opt_source((length(opt_source)-delay+1):length(opt_source)) = [];
        % Sum over all possible values that the two time-series may take.
        for i = unique(target_past')
            for j = unique(opt_source')
                % Record number of instances instead of probability.
                jointprob = sum(sum([opt_source target_past target_future]==[j i target_future_value],2)==3);
                % Discard cases with probability zero.
                if jointprob == 0
                    disp(['Pr(source,target_future,target_past=', num2str(j), ',', num2str(target_future_value), ',', num2str(i), ') is zero. Case discarded.'])
                else
                    futureprob = sum(target_future==target_future_value);
                    if futureprob == 0
                    else
                        pastprob = sum(sum([opt_source target_past]==[j i],2)==2);
                        if pastprob == 0
                        else
                        spec_info = spec_info + jointprob / futureprob * log(jointprob / pastprob / futureprob * length(target_future)) / log(2);
                        % Multiply term within log() by length to ensure dimensions are correct. Divide by log(2) to return units of bits.
                        end
                    end
                end
            end
        end
    % Optional time-series input is not given in this case.
    elseif nargin == 3
        % Create past and future time-series for target using given
        % time-delay by truncating at start or end of time-series.
        target_future = target;
        target_future(1:delay) = [];
        target_past = target;
        target_past((length(target)-delay+1):length(target)) = [];
        % Sum over all possible values that the time-series may take.
        for i = unique(target_past')
            jointprob = sum(sum([target_past target_future]==[i target_future_value],2)==2);
            % Discard cases with probability zero.
            if jointprob == 0
                disp(['Pr(target_past,target_future=', num2str(i), ',', num2str(target_future_value), ') is zero. Case discarded.'])
            else
                futureprob = sum(target_future==target_future_value);
                if futureprob == 0
                else
                    pastprob = sum(target_past==i);
                    if pastprob == 0
                    else
                        spec_info = spec_info + jointprob / futureprob * log(jointprob / pastprob / futureprob * length(target_future)) / log(2);
                        % Multiply term within log() by length to ensure dimensions are correct. Divide by log(2) to return units of bits.
                    end
                end
            end
        end
    else
        error('Number of arguments must be 3 or 4.')
    end
end