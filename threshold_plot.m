function fig_handle = threshold_plot(mat, mat25, mat50, mat75)

    fig_handle = figure;
    
    sub = subplot(2,2,1);
    histogram(mat(mat>0), 'Normalization', 'pdf');
    title('All recruitment weights for trial 0 test epoch');
    xlabel('bits');
    ylabel('pdf');
    sub.XLim = [-0.005 0.4];
    sub.YLim = [0 12];
    
    
    sub25 = subplot(2,2,2);
    histogram(mat25(mat25>0), 'Normalization', 'pdf');
    title('Thresholded at 25th percentile of functional weights');
    xlabel('bits');
    ylabel('pdf');
    sub25.XLim = [-0.005 0.4];
    sub25.YLim = [0 12];
    
    sub50 = subplot(2,2,3);
    histogram(mat50(mat50>0), 'Normalization', 'pdf');
    title('Thresholded at 50th percentile of functional weights');
    xlabel('bits');
    ylabel('pdf');
    sub50.XLim = [-0.005 0.4];
    sub50.YLim = [0 12];
    
    
    sub75 = subplot(2,2,4);
    histogram(mat75(mat75>0), 'Normalization', 'pdf');
    title('Thresholded at 75th percentile of functional weights');
    xlabel('bits');
    ylabel('pdf');
    sub75.XLim = [-0.005 0.4];
    sub75.YLim = [0 12];

end