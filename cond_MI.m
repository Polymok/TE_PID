% Given three time-series X, Y, Z, this function returns the conditional
% mutual information between X and Y conditioned on Z, denoted MI(X;Y|Z).
% Time-series should be given as column vectors.
%
% This function can only take discrete time-series.
%
% This function may take vector-valued time-series.

function [condMI] = cond_MI(X, Y, Z)
    % Ensure input time-series are column vectors.
    if size(X,1) < size(X,2)
        str = input('Input vector X has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            X = X';
        end
    end
    if size(Y,1) < size(Y,2)
        str = input('Input vector Y has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            Y = Y';
        end
    end
    if size(Z,1) < size(Z,2)
        str = input('Input vector Z has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            Z = Z';
        end
    end
    % Check if inputs are the same length.
    if (size(X,1)~=size(Y,1)) || (size(X,1)~=size(Z,1)) || (size(Y,1)~=size(Z,1))
        error('Input vectors are not of equal length.')
    end
    condMI = 0; % Initialize output.
    % Sum over all unique values of X, Y, and Z.
    for i = unique(X,'row')'
        for j= unique(Y,'row')'
            for k = unique(Z,'row')'
                % Record number of instances instead of probability.
                jointprob = sum(sum([X Y Z]==[i' j' k'],2)==(size(X,2)+size(Y,2)+size(Z,2)));
                % Discard cases with probability zero.
                if jointprob == 0
%                     disp(['Pr(X,Y,Z=', num2str(i'), ',', num2str(j'), ',', num2str(k'), ') is zero. Case discarded.'])
                else
                    probXZ = sum(sum([X Z]==[i' k'],2)==(size(X,2)+size(Z,2)));
                    probYZ = sum(sum([Y Z]==[j' k'],2)==(size(Y,2)+size(Z,2)));
                    if (probXZ == 0) || (probYZ == 0)
%                         disp(['Pr(X,Z=', num2str(i'), ',', num2str(k'), ') or Pr(Y,Z=', num2str(j'), ',', num2str(k'), ') is zero. Case discarded.'])
                    else
                        probZ = sum(Z==k');
                        condMI = condMI + jointprob/size(X,1) * log(jointprob * probZ / probXZ / probYZ) / log(2);
                    end    
                end
            end
        end
    end       
end