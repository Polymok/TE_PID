% Given three time-series X, Y, Z, this function returns the conditional
% mutual information between X and Y conditioned on Z, denoted MI(X;Y|Z).
%
% Output is a scalar in units of bits.
%
% This function can only take discrete, scalar-valued time-series.
%
% TODO:
% How to condition on events with probability zero?

function [condMI] = cond_MI(X, Y, Z)
    % Check if inputs are of acceptable type and length.
    if (~isvector(X)) || (~isvector(Y)) || (~isvector(Z))
        error('Inputs are not vectors.')
    elseif (length(X)~=length(Y)) || (length(X)~=length(Z)) || (length(Y)~=length(Z))
        error('Input vectors are not of equal length.')
    end
    % Ensure inputs are column vectors.
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    condMI = 0; % Initialize output.
    % The following three for loops sum over all possible values that the three time-series may take.
    for i = unique(X')
        for j = unique(Y')
            for k = unique(Z')
                jointprob = sum(sum([X Y Z]==[i j k],2)==3);
                % Cases of probability zero are discarded.
                if jointprob == 0
                    disp(['Pr(source,target_past,target_future=', num2str(i), ',', num2str(k), ',', num2str(j), ') is zero. Case discarded.'])
                else
                    prob1condition = sum(sum([X Z]==[i k],2)==2);
                    prob2condition = sum(sum([Y Z]==[j k],2)==2);
                    if (prob1condition == 0) || (prob2condition == 0)
                        disp(['Pr(source,target_past=', num2str(i), ',', num2str(k), ') or Pr(target_past,target_future=', num2str(k), ',', num2str(j), ') is zero. Case discarded.'])
                    else
                        probcondition = sum(Z==k);
                        condMI = condMI + jointprob/length(X) * log(jointprob * probcondition / prob1condition / prob2condition) / log(2);
                    end    
                end
            end
        end
    end
end