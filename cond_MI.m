% Given three time-series X, Y, Z, this function returns the conditional
% mutual information between X and Y conditioned on Z, denoted MI(X;Y|Z).
% Time-series should be given as column vectors. Output is given in units
% of bits.
%
% This function can only take discrete time-series.

function [condMI] = cond_MI(X, Y, Z)

%     %% Check inputs.
%     % Ensure input time-series are column vectors.
%     if size(X,1) < size(X,2)
%         str = input('Input vector X has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
%         if str == 'y'
%             X = X';
%         end
%     end
%     if size(Y,1) < size(Y,2)
%         str = input('Input vector Y has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
%         if str == 'y'
%             Y = Y';
%         end
%     end
%     if size(Z,1) < size(Z,2)
%         str = input('Input vector Z has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
%         if str == 'y'
%             Z = Z';
%         end
%     end
%     if (size(X,1)~=size(Y,1)) || (size(X,1)~=size(Z,1)) || (size(Y,1)~=size(Z,1))
%         error('Input vectors are not of equal length.')
%     end
    
    %% Calculate condition mutual information.
    condMI = 0; % Initialize output.
    % Sum over all unique values of X, Y, and Z.
    for i = unique(X,'row')' % Specifying 'row' argument is necessary in case inputs are vectors. Taking the transpose here is necessary to facilitate the for loop. Re-take transpose within for loop.
        for j= unique(Y,'row')'
            for k = unique(Z,'row')'
                % Count number of instances instead of probability.
                probXYZ = sum(sum([X Y Z]==[i' j' k'],2)==(size(X,2)+size(Y,2)+size(Z,2))); % p(x,y,z)
                probXZ = sum(sum([X Z]==[i' k'],2)==(size(X,2)+size(Z,2))); % p(x,z)
                probYZ = sum(sum([Y Z]==[j' k'],2)==(size(Y,2)+size(Z,2))); % p(y,z)
                % Discard cases with probability zero and print message to command window.
                if probXYZ == 0
%                     disp(['Pr(X,Y,Z=', num2str(i'), ',', num2str(j'), ',', num2str(k'), ') is zero. Case discarded.'])
                elseif probXZ == 0
%                     disp(['Pr(X,Z=', num2str(i'), ',', num2str(k'), ') is zero. Case discarded.'])
                elseif probYZ == 0
%                     disp(['Pr(Y,Z=', num2str(j'), ',', num2str(k'), ') is zero. Case discarded.'])
                else
                    probZ = sum(Z==k'); % p(z)
                    condMI = condMI + probXYZ/size(X,1) * log2(probXYZ * probZ / probXZ / probYZ); % MI(X;Y|Z) = sum_{X=x,Y=y,Z=z} p(x,y,z)*log(p(x,y,z)*p(z)/p(x,z)/p(y,z))
                end
            end
        end
    end   
    
end