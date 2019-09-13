% Plot histograms for partial information terms from all triplets,
% functional network, and recruitment network.
%
% Inputs must be nx7 matrices, e.g. outputs of TE_PID.m.

function fig_handle =  PID_plot(varargin)
    fig_handle = figure;
    for i = 1:nargin
        subplot(nargin,3,(i-1)*3+1);
        histogram(varargin{i}(:,4), 'Normalization', 'pdf');
        title('Synergy');
    
        subplot(nargin,3,(i-1)*3+2);
        histogram(varargin{i}(:,5), 'Normalization', 'pdf');
        title('Redundancy');

        subplot(nargin,3,(i-1)*3+3);
        histogram(varargin{i}(:,6:7), 'Normalization', 'pdf');
        title('Unique');
    end
end