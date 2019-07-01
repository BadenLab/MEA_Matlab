function Kernels = calculate_Kernels_speed (Epochs, Input, kbinsize, time_window)

    spiketimestamps = gather(Input);
    stx = length(spiketimestamps(1,:));
    sty = length(spiketimestamps(:,1));
    post_window = 0.2;
    
    time_window_temp = time_window;
    time_window = time_window+post_window;
    time_window_in_ms = time_window*1000;
    kernel_bins = (time_window/kbinsize);
    kernel_bins_idx = int16(kernel_bins);
    %Check which noise stimuli are available:
    
    N_stimcode(1) = isfield(Epochs, 'Noise380');
    N_stimcode(2) = isfield(Epochs, 'Noise430');
    N_stimcode(3) = isfield(Epochs, 'Noise480');
    N_stimcode(4) = isfield(Epochs, 'Noise505');
    N_stimcode(5) = isfield(Epochs, 'Noise560');
    N_stimcode(6) = isfield(Epochs, 'Noise630');
    
    Noise_size = nnz(N_stimcode);
    
    
    
    Kernels = zeros(Noise_size,kernel_bins_idx,stx);
    
    Noise_time = Epochs.Epochs_duration_o(1,Epochs.CNoise)/6000;
    Noise_freq = 1/Noise_time;
%      Noise_freq = str2double(inputdlg(prompt,title,dims,definput));
    time_per_frame = 1/Noise_freq;
     
     
    C_epochs = Epochs.CNoise;
     %Calculate the number of colour noise epochs
    nr_C_epochs = nnz(C_epochs);
     %Check how long the colour noise stimulus is running in s
    sequence_time = unique(ceil(Epochs.Epochs_end(C_epochs)-Epochs.Epochs_begin(C_epochs)));
%      sequence_time = sequence_time;
     Noise_seq_complete = NaN(6000,Noise_size);
     Repeats_per_epoch = 1;
     
     
     %get the stimulus data as logicals and according to the number of
     %repeats
     nidx = 1;

     for nn = 1:length(N_stimcode)
         
     if N_stimcode(nn) == 1
         switch nn
             case 1
                 log_Noise380 = logical(Epochs.Noise380);
                 log_Noise380 = repmat(log_Noise380,Repeats_per_epoch);
                 log_Noise380 = log_Noise380(:,1);
                 Noise_seq_complete(:,nidx) = log_Noise380;
                 nidx = nidx+1;
                 
             case 2
                 log_Noise430 = logical(Epochs.Noise430);
                 log_Noise430 = repmat(log_Noise430,Repeats_per_epoch);
                 log_Noise430 = log_Noise430(:,1);
                 Noise_seq_complete(:,nidx) = log_Noise430;
                 nidx = nidx+1;
             case 3
                 log_Noise480 = logical(Epochs.Noise480);
                 log_Noise480 = repmat(log_Noise480,Repeats_per_epoch);
                 log_Noise480 = log_Noise480(:,1);
                 Noise_seq_complete(:,nidx) = log_Noise480;
                 nidx = nidx+1;
             case 4
                 log_Noise505 = logical(Epochs.Noise505);
                 log_Noise505 = repmat(log_Noise505,Repeats_per_epoch);
                 log_Noise505 = log_Noise505(:,1);
                 Noise_seq_complete(:,nidx) = log_Noise505;
                 nidx = nidx+1;
             case 5
                 log_Noise560 = logical(Epochs.Noise560);
                 log_Noise560 = repmat(log_Noise560,Repeats_per_epoch);
                 log_Noise560 = log_Noise560(:,1);
                 Noise_seq_complete(:,nidx) = log_Noise560;
                 nidx = nidx+1;
             case 6
                 log_Noise630 = logical(Epochs.Noise630);
                 log_Noise630 = repmat(log_Noise630,Repeats_per_epoch);
                 log_Noise630 = log_Noise630(:,1);
                 Noise_seq_complete(:,nidx) = log_Noise630;
                 nidx = nidx+1;
         end
     end
     
     end
     
  Upsample_factor = (1/Noise_freq)/kbinsize;
a = 1;
for ii = 1:length(Noise_seq_complete(:,1))
    for tt = 1:Upsample_factor
        Noise_seq_upsample(a,:) = Noise_seq_complete(ii,:);
        a = a+1;
    end
end
% Add time window before and after the set time window
pre_time_ms = time_window_temp*1000;
post_time_ms = post_window * 1000;
idx_end = length(Noise_seq_upsample(:,1));
Noise_seq_upsample_temp(1:pre_time_ms,4) = zeros;
Noise_seq_upsample_temp(pre_time_ms+1:idx_end+pre_time_ms,:) = Noise_seq_upsample;
idx_end = length(Noise_seq_upsample_temp(:,1));
Noise_seq_upsample_temp(idx_end+1:idx_end+post_time_ms,:) = 0;
Noise_seq_upsample = Noise_seq_upsample_temp;
clear Noise_seq_upsample_temp


% Create a vector with the timings for the stimulus miliseconds


     
     colour_begin = Epochs.Epochs_begin(C_epochs);
     colour_end = Epochs.Epochs_end(C_epochs);
     true_spikes_b_e = zeros(2,stx);
     stim_ms = (colour_begin:0.001:colour_end+0.001);
     
     %% Loop
     
parfor cc = 1:stx
%          ff = cc/stx;
%             waitbar(ff);
     cc
     %Check how many spikes are in the cluster
     temp_spiketimestamps = spiketimestamps(:,cc);
     %temp_spiketimestamps(temp_spiketimestamps == 0) = NaN;
     %nan_spiketimestamps = isnan(temp_spiketimestamps);
     %l_spiketimestamps = length(nan_spiketimestamps) - nnz(nan_spiketimestamps);
     
     %Try to find the position of the first spike inside the colour noise
   
     Begin = 1;
     End = sty;
     test_begin = 0;
     while test_begin == 0
     test_begin = temp_spiketimestamps(Begin,1) > colour_begin;
     if temp_spiketimestamps(Begin,1) == 0
             break
     end
     Begin = Begin + 100;
     
     end
     test_begin = 0;
     while test_begin == 0
         Begin = Begin -1;
         if Begin == 0
            break
         end
         test_begin = temp_spiketimestamps(Begin,1)< colour_begin;
         if isnan(temp_spiketimestamps(Begin,1))
             break
         end
     end
     if Begin == 0
         continue
     else
     Begin = Begin +1;
     end
     % Find the end of the true spikes
     test_end = 0;
     while test_end == 0
         if temp_spiketimestamps(End,1) == 0
             temp_spiketimestamps(End,1) = NaN;
         end
     test_end = temp_spiketimestamps(End,1) < colour_end;
     
     End = End - 100;
     if End <= 1
         break
     end
     
     end
     test_end = 0;
     while test_end == 0
         End = End +1;
     if temp_spiketimestamps(End,1) == 0
          break
     end
         test_end = temp_spiketimestamps(End,1) > colour_end;
     if End <= 1
         break
     end
     
     end
     End = End-1;
     
     true_spikes_b_e(:,cc) = [Begin End]';
     
     
end

%Check which traces are empty
empty_true_spikes = find(~true_spikes_b_e(1,:));

stx_new = (1:1:stx);
stx_new(empty_true_spikes) = [];
     
%% Loop 2

parfor kk = 1:stx
    skip = empty_true_spikes;
    if kk == skip
        continue
    else
    spiketimestamps_temp = spiketimestamps(:,kk);
    Begin_End_real = true_spikes_b_e(:,kk);
    real_spikes = spiketimestamps_temp(Begin_End_real(1,1):Begin_End_real(2,1),1);
    stim_time = stim_ms;
    spike_history = zeros(length(real_spikes),kernel_bins,4);
    idx_spike = 1;
    for ii = 1:length(real_spikes(:,1))
        spike = real_spikes(ii,1);
        test_spike = 0;
       
        while test_spike == 0
            test_spike = spike<stim_time(idx_spike);
            idx_spike = idx_spike+1;
        end
        idx_spike = idx_spike-2;
        Noise_seq_upsample_temp = Noise_seq_upsample;
        spike_history(ii,:,:) = Noise_seq_upsample_temp(idx_spike:idx_spike+time_window_in_ms-1,:);
        
    end
    Kernels(:,:,kk) = squeeze(sum(spike_history,1))';
    kk
    end
    end
end
        
            
        
        


    
     
    
     
     
     
     
     
     
     
     
%      
%      i1 = 1;
%      i2 = 0;
%      clear true_spikes
%      for ee = 1:nr_C_epochs
%          
%          
%          true_spikes_temp = (temp_spiketimestamps >= colour_begin(ee)) .* ...
%              (temp_spiketimestamps <= colour_end(ee));
%          true_spikes_temp = logical(true_spikes_temp);
%          %count how many spikes were found
%          i2 = i2 + nnz(true_spikes_temp);
%         true_spikes(i1:i2,1) = temp_spiketimestamps(true_spikes_temp);
%         true_spikes(i1:i2) = true_spikes(i1:i2) - colour_begin(ee);
%         i1 = length(true_spikes) + 1;
%         %Set the beginning of the noise stimulus as spiketime 0
%         
%      end
   
      
   %length of new spiketimestamps matrix 
%      end
% 
%  tsy = length(true_spikes(:,1));


