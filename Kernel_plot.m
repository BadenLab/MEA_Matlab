for ii = KQI_past% find(Best_cells')%KQI_past %5%1:100%100:200
    Kernel_plot = squeeze(gather(Epochs.Kernels(:,:,ii)));
    
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
     figure
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
    hold off
end
    


