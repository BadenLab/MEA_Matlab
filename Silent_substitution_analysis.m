for ii = [7,8,9,10]
    SS_max(ii,:) = max(gather(Find_Epoch(Epochs,ii)));
end

[~,SS_max_max] = max(SS_max); 

SS_response(1,:) = SS_max_max == 7;
SS_response(2,:) = SS_max_max == 8;
SS_response(3,:) = SS_max_max == 9;
SS_response(4,:) = SS_max_max == 10;


colourstring = [1 0 1; 0 0 1; 0 1 0; 1 0 0];
for kk = 1:4
    centres_temp = centres(SS_response(kk,:),:);
    
    scatter(centres_temp(:,1),centres_temp(:,2),'MarkerFaceColor',colourstring(kk,:), 'MarkerEdgeColor', 'None')
    hold on
    
end

hold off



for kk = 4
    
    plot_temp = normalize(mean(Smoothened_averaged_spikes(:,SS_response(kk,:)),2));
    
    plot(plot_temp,'Color',colourstring(kk,:))
    
    hold on
    
    plot_temp1 = Smoothened_averaged_spikes(:,SS_response(kk,:));
    
    plot_temp2 = plot_temp + plot_temp1;
    plot_temp3 = plot_temp - plot_temp1;
    
    
   
    a = area(plot_temp2);
    a.FaceColor = [0.9 0.9 0.9];
    b = area(plot_temp3,'k');
    b.FaceColor = [0.9 0.9 0.9];
    alpha(a,.5)
    alpha(b,.5)
        
        
    
end

