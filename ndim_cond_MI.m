% Given three time-series A, B, C, this function returns the conditional
% mutual information between A and B conditioned on C, denoted MI(A;B|C).
%
% This function can only take discrete time-series.
%
% This function may take vector-valued time-series.
%
% STATUS: incomplete.

function [condMI] = ndim_cond_MI(series1, series2, condition)
    if (length(series1)~=length(series2)) || (length(series1)~=length(condition)) || (length(series2)~=length(condition))
        error('Input vectors are not of equal length.')
    end
    if size(series1,1) < size(series1,2)
        str = input('Input vector A has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            series1 = series1';
        end
    end
    if size(series2,1) < size(series2,2)
        str = input('Input vector B has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            series2 = series2';
        end
    end
    if size(condition,1) < size(condition,2)
        str = input('Input vector C has a greater number of columns than rows. Each column should contain the entire time-series of a single neuron. Transpose input matrix? y/n: ','s');
        if str == 'y'
            condition = condition';
        end
    end
    timematrix = [series1 series2 condition];
    condMI = 0;
    for I = 1:size(series1,2)
        for J = 1:size(series2,2)
            for K = 1:size(condition,2)
                for i = unique(series1)
                    for j = unique(series2)
                        for k = unique(condition)
                            jointprob = sum(sum(timematrix==[i j k],2)==size(timematrix,2))/size(timematrix,1);
                            prob2condition = sum(sum([timematrix(:,2) timematrix(:,3)]==[j k],2)==2);
                            if prob2condition == 0
                                prob1_2condition = 0;
                            else
                                 prob1_2condition = sum(sum(timematrix==[i j k],2)==size(timematrix,2))/prob2condition;
                            end
                            prob1_condition = sum(sum([timematrix(:,1) timematrix(:,3)]==[i k],2)==2)/sum(timematrix(:,3)==k);
                            condMI = condMI + jointprob * log(prob1_2condition / prob1_condition) / log(2);
                        end
                    end
                end
            end
        end
    end
end