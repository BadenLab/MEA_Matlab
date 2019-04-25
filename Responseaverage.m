function [Epochs_average, Epochs] = Responseaverage (Epochs, Bined_epochs, Mode)
if Mode == 0
for ii = 1: Epochs.nr_epochs
    
    stimuli_nr = length(Epochs.stimulus_starts(:,1,1));
    
    
    Epochs_bined_temp = squeeze(Bined_epochs(:,:,:,ii));
    Epochs_average_temp = squeeze(nanmean(Epochs_bined_temp,1));
    Epochs_average(:,:,ii) = Epochs_average_temp;
end
else
   
    %identify the unique epochs
    
    %Check for epochs that were not repeated but shown only once
    
    unique_ep = unique(Epochs.Epoch_code);
    l_unique_ep = length(unique_ep);
    nr_unique = NaN(1,l_unique_ep);
    
    for ss = 1:l_unique_ep
       nr_unique(1,ss) = sum(Epochs.Epoch_code(1,:) == unique_ep(ss));
        
    end
    
    epochs_one_repeat = nr_unique == 1;
    nr_epochs_one_repeat = nnz(epochs_one_repeat);

    
    e_nor = epochs_one_repeat == 0;
    
    unique_ep_n1 = unique_ep(e_nor);
    l_unique_ep_n1 = length(unique_ep_n1);
    
    nr_ep = max(nr_unique);
     
    ep_idx = NaN(l_unique_ep_n1,nr_ep);
    
    Colour_noise_epochs = Epochs.CNoise == 0;
    Colour_noise_epochs = Colour_noise_epochs(1,1:l_unique_ep_n1);
    Colour_noise_epochs = find(Colour_noise_epochs);
    Bined_epochs = squeeze(Bined_epochs(1,:,:,:));
    
    
    
    
    
    for ii = 1:l_unique_ep_n1
        aa = unique_ep_n1(ii);
        ep_idx(ii,:) = find(Epochs.Epoch_code == aa);
    end
    %sort the epochs so that they are in the sequence in which 
    %they have been shown
    [~,idx] = sort(ep_idx(:,1));
    ep_idx = ep_idx(idx,:);
    l_B_x = length(Bined_epochs(:,1,1));
    l_B_y = length(Bined_epochs(1,:,1));
    
    

    mean_bined_epochs = NaN(l_B_y,length(ep_idx(:,1)),l_B_x);
    
    for ii = 1:l_unique_ep_n1
        for kk = 1:length(Bined_epochs(:,1,1))
        temp_Bined_epochs = squeeze(Bined_epochs(kk,:,ep_idx(ii,:)));
        mean_bined_epochs(:,ii,kk) = nanmean(temp_Bined_epochs,2);
        end
    end
%         test_mean = mean_bined_epochs(:,:,89);
        
%         test_mean_whole = NaN(1,(max(Epochs.nr_bins(1,1,:))));
        aa = 1;
        
    for ii = Colour_noise_epochs
         bb = Epochs.nr_bins(1,1,ii)+(aa-1);
        for kk = 1:length(mean_bined_epochs(1,1,:))
       
        mean_bined_whole((aa:bb),kk) = mean_bined_epochs((1:Epochs.nr_bins(1,1,ii)),ii,kk);
        
        end
        aa = bb+1;
    end
    
    Epochs_average = mean_bined_whole;
    Epochs.loop_repeats = nr_ep;
    
        
        
        
        
        
        
        
        
        
    
    
    
    
end
    
    
    
    





end