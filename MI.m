% Given two vectors X and Y, this function returns the mutual information
% between X and Y.
%
% Output is a scalar in units of bits.
%
% This function can only take discrete, scalar-valued vectors.

function [mutual_info] = MI(X, Y)
    % Check if inputs are of acceptable type and length.
    if (~isvector(X)) || (~isvector(Y))
        error('Inputs must be vectors.')
    elseif (length(X)~=length(Y))
        error('Input vectors are not of equal length.')
    end
    % Ensure inputs are column vectors.
    X = X(:);
    Y = Y(:);
    mutual_info = 0; % Initialize output.
    % Sum over all possible values that the two time-series may take.
    for i = unique(X')
        for j = unique(Y')
            % Record number of instances instead of probability.
            jointprob = sum(sum([X Y]==[i j],2)==2);
            % Discard cases with probability zero.
            if jointprob == 0
                disp(['Pr(X,Y=', num2str(i), ',', num2str(k), ') is zero. Case discarded.'])
            else
                probx = sum(X==i);
                proby = sum(Y==j);
                mutual_info = mutual_info + jointprob / length(X) * log(jointprob / probx / proby * length(X)) / log(2);
                % Divide number of instances by length to obtain probability. Divide by log(2) to return units of bits.
            end
        end
    end
end