%CLuster Plotting
fig = figure
fig.Renderer='Painters';
ax1 = subplot(11,3,1);
stdshade(Traces_Cluster_1',0.5,'k');
set(ax1,'XTick',[])

ax2 = subplot(11,3,4);
stdshade(Traces_Cluster_2',0.5,'k');
set(ax2,'XTick',[])

ax3 = subplot(11,3,7);
plot(Traces_Cluster_3,'k');
set(ax3,'XTick',[])

ax4 = subplot(11,3,10);
stdshade(Traces_Cluster_4',0.5,'k');
set(ax4,'XTick',[])

ax5 = subplot(11,3,13);
stdshade(Traces_Cluster_5',0.5,'k');
set(ax5,'XTick',[])

ax6 = subplot(11,3,16);
stdshade(Traces_Cluster_6',0.5,'k');
set(ax6,'XTick',[])

ax7 = subplot(11,3,19);
plot(Traces_Cluster_7,'k');
set(ax7,'XTick',[])

ax8 = subplot(11,3,22);
stdshade(Traces_Cluster_8',0.5,'k');
set(ax8,'XTick',[])

ax9 = subplot(11,3,25);
stdshade(Traces_Cluster_9',0.5,'k');
set(ax9,'XTick',[])

ax10 = subplot(11,3,28);
stdshade(Traces_Cluster_10',0.5,'k');
set(ax10,'XTick',[])

ax11 = subplot(11,3,31);
stimulus_plot = bar(temp_sfg,2,'FaceColor','flat','EdgeColor','flat');
        for kk = 1:length(RGB_whole_trace)
    
        stimulus_plot.CData(kk,:) = RGB_whole_trace(kk,:);
        end
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11],'x')

ax12 = subplot(11,3,2);
h1 = heatmap(Traces_Cluster_1','GridVisible','off','Colormap',c);
h1.ColorbarVisible = 'off';

ax13 = subplot(11,3,5);
h2 = heatmap(Traces_Cluster_2','GridVisible','off','Colormap',c);
h2.ColorbarVisible = 'off';

ax14 = subplot(11,3,8);
h3 = heatmap(Traces_Cluster_3','GridVisible','off','Colormap',c);
h3.ColorbarVisible = 'off';

ax15 = subplot(11,3,11);
h4 = heatmap(Traces_Cluster_4','GridVisible','off','Colormap',c);
h4.ColorbarVisible = 'off';

ax16 = subplot(11,3,14);
h5 = heatmap(Traces_Cluster_5','GridVisible','off','Colormap',c);
h5.ColorbarVisible = 'off';

ax17 = subplot(11,3,17);
h4 = heatmap(Traces_Cluster_6','GridVisible','off','Colormap',c);
h4.ColorbarVisible = 'off';

ax18 = subplot(11,3,20);
h5 = heatmap(Traces_Cluster_7','GridVisible','off','Colormap',c);
h5.ColorbarVisible = 'off';

ax19 = subplot(11,3,23);
h6 = heatmap(Traces_Cluster_8','GridVisible','off','Colormap',c);
h6.ColorbarVisible = 'off';

ax20 = subplot(11,3,26);
h7 = heatmap(Traces_Cluster_9','GridVisible','off','Colormap',c);
h7.ColorbarVisible = 'off';

ax21 = subplot(11,3,29);
h8 = heatmap(Traces_Cluster_10','GridVisible','off','Colormap',c);
h8.ColorbarVisible = 'off';

ax22 = subplot(11,3,32);
stimulus_plot = bar(temp_sfg,2,'FaceColor','flat','EdgeColor','flat');
        for kk = 1:length(RGB_whole_trace)
    
        stimulus_plot.CData(kk,:) = RGB_whole_trace(kk,:);
        end
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11],'x')

ax23 = subplot(11,3,3);
    Kernel_plot = squeeze(Cluster_mean(:,:,1));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
ax24 = subplot(11,3,6);
    Kernel_plot = squeeze(Cluster_mean(:,:,2));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
ax25 = subplot(11,3,9);
    Kernel_plot = squeeze(Cluster_mean(:,:,3));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
   
   ax26 = subplot(11,3,12);
    Kernel_plot = squeeze(Cluster_mean(:,:,4));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
   ax27 = subplot(11,3,15);
    Kernel_plot = squeeze(Cluster_mean(:,:,5));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
   ax28 = subplot(11,3,18);
    Kernel_plot = squeeze(Cluster_mean(:,:,6));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
   
   ax29 = subplot(11,3,21);
    Kernel_plot = squeeze(Cluster_mean(:,:,7));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off

   
   ax30 = subplot(11,3,24);
    Kernel_plot = squeeze(Cluster_mean(:,:,8));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
   
   ax31 = subplot(11,3,27);
    Kernel_plot = squeeze(Cluster_mean(:,:,9));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
   ax35 = subplot(11,3,30);
    Kernel_plot = squeeze(Cluster_mean(:,:,10));
    colorstring = 'mbgr';
    xaxis_k = (-1.999:0.001:0.2);
    for kk = 1:length(Kernel_plot(:,1))
    Kernel_plot = Kernel_plot./max(max(Kernel_plot));
    plot(xaxis_k,Kernel_plot(kk,:),'Color',colorstring(kk),'LineWidth',2)
    
    hold on
    end
%     line([0 0], [-5 5],'Color','k'); 
    xlim ([-0.5 0.2])
%     set(gca, 'YTick', [])
   hold off
   
   
