% Wrapper function to find intersection between list of functional or
% recruitment triplets and the list of all triplets, as well as
% corresponding PID terms.

function subset_PID = PID_extract(all_PID, list)

%     if size(all_PID,2)~=7
%         error('Input matrix of all PID values must have 7 columns.')
%     elseif size(list,2)~=3
%         error('Input list of triplets must have 3 columns.')
%     end
    [~,iAll] = intersect(all_PID(:,1:3), list, 'rows');
    subset_PID = all_PID(iAll,:);
    
end
