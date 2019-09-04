% Plot histograms for partial information terms from all triplets,
% functional network, and recruitment network.

function [fig_handle, synergy, redundancy, unique] =  PID_plot(all, functional, recruitment)
    if nargin ~= 3
        error('3 inputs required.')
    elseif ~ismatrix(all) || ~ismatrix(functional) || ~ismatrix(recruitment)
        error('Inputs must be matrices.')
    elseif size(all,2)~=7 || size(functional,2)~=7 || size(recruitment,2)~=7
        error('Input matrices must have 7 columns.')
    end
    
    fig_handle = figure;
    
    subplot(3,3,1);
    synergy_all = histogram(all(:,4), 'Normalization', 'pdf');
    title('Synergy: all triplets');
    
    subplot(3,3,4);
    synergy_func = histogram(functional(:,4), 'Normalization', 'pdf');
    title('Synergy: functional network');
    
    subplot(3,3,7);
    synergy_recr = histogram(recruitment(:,4), 'Normalization', 'pdf');
    title('Synergy: recruitment network');
    
    subplot(3,3,2);
    redundancy_all = histogram(all(:,5), 'Normalization', 'pdf');
    title('Redundancy: all triplets');
    
    subplot(3,3,5);
    redundancy_func = histogram(functional(:,5), 'Normalization', 'pdf');
    title('Redundancy: functional network');
    
    subplot(3,3,8);
    redundancy_recr = histogram(recruitment(:,5), 'Normalization', 'pdf');
    title('Redundancy: recruitment network');
    
    subplot(3,3,3);
    unique_all = histogram(all(:,6:7), 'Normalization', 'pdf');
    title('Unique: all triplets');
    
    subplot(3,3,6);
    unique_func = histogram(functional(:,6:7), 'Normalization', 'pdf');
    title('Unique: functional network');
    
    subplot(3,3,9);
    unique_recr = histogram(recruitment(:,6:7), 'Normalization', 'pdf');
    title('Unique: recruitment network');
    
    synergy = [synergy_all synergy_func synergy_recr];
    redundancy = [redundancy_all redundancy_func redundancy_recr];
    unique = [unique_all unique_func unique_recr];
end