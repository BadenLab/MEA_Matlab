
for ii = 1:length(KQI_past(1,:))
    Kernel_UV(ii,:) = squeeze(Kernels(1,:,KQI_past(1,ii)));
    Kernel_Blue(ii,:) = squeeze(Kernels(2,:,KQI_past(1,ii)));
    Kernel_Green(ii,:) = squeeze(Kernels(3,:,KQI_past(1,ii)));
    Kernel_Red(ii,:) = squeeze(Kernels(4,:,KQI_past(1,ii)));
end

for ii = 1:273
    Kernel_UV(ii,:) = squeeze(Kernels(1,:,ii));
    Kernel_Blue(ii,:) = squeeze(Kernels(2,:,ii));
    Kernel_Green(ii,:) = squeeze(Kernels(3,:,ii));
    Kernel_Red(ii,:) = squeeze(Kernels(4,:,ii));
end
    [Cluster_UV,C_UV] = kmeans(Kernel_UV,5,'Distance','correlation',...
    'Replicates',500);
    [Cluster_Blue,C_Blue] = kmeans(Kernel_Blue,5,'Distance','correlation',...
    'Replicates',500);
    [Cluster_Green,C_Green] = kmeans(Kernel_Green,5,'Distance','correlation',...
    'Replicates',500);
    [Cluster_Red,C_Red] = kmeans(Kernel_Red,5,'Distance','correlation',...
    'Replicates',500);   
    



Cluster_master(:,1) = Cluster_UV;
Cluster_master(:,2) = Cluster_Blue;
Cluster_master(:,3) = Cluster_Green;
Cluster_master(:,4) = Cluster_Red;

[Cluster_m_idx, C_m_idx] = kmeans(Cluster_master,5,'Replicates',500);   

Cluster_m = nan(4,2200,length(KQI_past(1,:)));
Cluster_m = nan(4,2200,273);
Cluster_mean = nan(4,2200,3);
for ii = 1:5
    Cluster_log = Cluster_m_idx == ii;
    nnz_Cluster_log = nnz(Cluster_log);
    for kk = 1:4
        Cluster_m(kk,:,1:nnz_Cluster_log) = Kernels(kk,:,KQI_past(1,Cluster_log));
        Cluster_mean(kk,:,ii) = nanmean(Cluster_m(kk,:,:),3);
    end
end
 
for ii = 1:5
    Cluster_log = Cluster_m_idx == ii;
    nnz_Cluster_log = nnz(Cluster_log);
    for kk = 1:4
        Cluster_m(kk,:,1:nnz_Cluster_log) = Kernels(kk,:,Cluster_log);
        Cluster_mean(kk,:,ii) = nanmean(Cluster_m(kk,:,:),3);
    end
end
    

Traces_Cluster_1 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 1));
Traces_Cluster_2 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 2));
Traces_Cluster_3 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 3));
Traces_Cluster_4 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 4));
Traces_Cluster_5 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 5));
Traces_Cluster_6 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 6));
Traces_Cluster_7 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 7));
Traces_Cluster_8 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 8));
Traces_Cluster_9 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 9));
Traces_Cluster_10 = Smoothened_averaged_spikes(:,KQI_past(Cluster_m_idx == 10));

Traces_cluster_mean(:,1) = mean(Traces_Cluster_1,2);
Traces_cluster_mean(:,2) = mean(Traces_Cluster_2,2);
Traces_cluster_mean(:,3) = mean(Traces_Cluster_3,2);
Traces_cluster_mean(:,4) = mean(Traces_Cluster_4,2);
Traces_cluster_mean(:,5) = mean(Traces_Cluster_5,2);
Traces_cluster_mean(:,6) = mean(Traces_Cluster_6,2);
Traces_cluster_mean(:,7) = mean(Traces_Cluster_7,2);
Traces_cluster_mean(:,8) = mean(Traces_Cluster_8,2);
Traces_cluster_mean(:,9) = mean(Traces_Cluster_9,2);
Traces_cluster_mean(:,10) = mean(Traces_Cluster_10,2);

%%
Traces_Cluster_1 = Smoothened_averaged_spikes(:,Cluster_m_idx == 1);
Traces_Cluster_2 = Smoothened_averaged_spikes(:,Cluster_m_idx == 2);
Traces_Cluster_3 = Smoothened_averaged_spikes(:,Cluster_m_idx == 3);
Traces_Cluster_4 = Smoothened_averaged_spikes(:,Cluster_m_idx == 4);
Traces_Cluster_5 = Smoothened_averaged_spikes(:,Cluster_m_idx == 5);


Traces_cluster_mean(:,1) = mean(Traces_Cluster_1,2);
Traces_cluster_mean(:,2) = mean(Traces_Cluster_2,2);
Traces_cluster_mean(:,3) = mean(Traces_Cluster_3,2);
Traces_cluster_mean(:,4) = mean(Traces_Cluster_4,2);
Traces_cluster_mean(:,5) = mean(Traces_Cluster_5,2);
