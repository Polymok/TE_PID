% Given one 2-dimensional time-series X and two 1-dimensional time-series
% Y, Z, this function returns the conditional mutual information between X
% and Y conditioned on Z, denoted MI(X;Y|Z).
%
% Output is a scalar in units of bits.
%
% This function can only take discrete time-series. The first input
% time-series must be 2-dimensional. The second and third must be
% 1-dimensional. This design is particular to calculating the transfer
% entropy partial information decomposition.
%
% TODO:
% How to condition on events with probability zero?

function [condMI] = cond_MI_2dim(X, Y, Z)
    % Check if inputs are of acceptable type and length.
    if (~isvector(Y)) || (~isvector(Z))
        error('Y and Z are not vectors.')
    elseif (size(X,1)~=2) && (size(X,2)~=2)
        error('X is not 2-dimensional.')
    elseif (length(X)~=length(Y)) || (length(X)~=length(Z)) || (length(Y)~=length(Z))
        error('Input time-series are not of equal length.')
    end
    % Ensure that time-series are represented as column vectors.
    if size(X,1) < size(X,2)
        str = input('Input 2-dimensional time-series has greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input time-series? y/n: ','s');
        if str == 'y'
            X = X';
        end
        clear str
    end
    % Ensure inputs are column vectors.
    Y = Y(:);
    Z = Z(:);
    condMI = 0;
    % The following four for loops sum over all possible values that the
    % three time-series---one of which is 2-dimensional---may take.
    for i = unique(X(:,1)')
        for j = unique(X(:,2)')
            for k = unique(Y')
                for l = unique(Z')
                    jointprob = sum(sum([X(:,1) X(:,2) Y Z]==[i j k l],2)==4);
                    % Cases of probability zero are discarded.
                    if jointprob == 0
                        disp(['Pr(X_1,X_2,Y,Z=', num2str(i), ',', num2str(j), ',', num2str(k), ',', num2str(l), ') is zero. Case discarded.'])
                    else
                        probXZ = sum(sum([X(:,1) X(:,2) Z]==[i j l],2)==3);
                        probYZ = sum(sum([Y Z]==[k l],2)==2);
                        if (probXZ == 0) || (probYZ == 0)
                            disp(['Pr(X_1,X_2,Z=', num2str(i), ',', num2str(j), ',', num2str(l), ') or Pr(Y,Z=', num2str(k), ',', num2str(l), ') is zero. Case discarded.'])
                        else
                            probZ = sum(Z==l);
                            condMI = condMI + jointprob/size(X,1) * log(jointprob * probZ / probXZ / probYZ) / log(2);
                        end
                    end
                end
            end
        end
    end
end