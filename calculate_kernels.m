function Kernels = calculate_kernels(Epochs, kbinsize, time_window)
    spiketimestamps = gather(Epochs.spiketimestamps);
    stx = length(spiketimestamps(1,:));
    sty = length(spiketimestamps(:,1));
    post_window = 0.2;
    time_window = time_window+post_window;
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
    
    
    % Ask which frequency was used for the noise
    prompt = {"Enter the Frequency of the noise in Hz"};
     title = 'Frequency';
     dims = [1 35];
     definput = {'Frequency in Hz'};
     Noise_freq = str2double(inputdlg(prompt,title,dims,definput));
     time_per_frame = 1/Noise_freq;
     
     
     %% Find the epochs that are colour noise
     
     C_epochs = Epochs.Epoch_code == single(0.35);
     %Calculate the number of colour noise epochs
     nr_C_epochs = nnz(C_epochs);
     %Check how long the colour noise stimulus is running
     sequence_time = unique(ceil(Epochs.Epochs_end(C_epochs)-Epochs.Epochs_begin(C_epochs)));
%      sequence_time = sequence_time;
     
% Output error when the time of the colour noise stimulus is not consistent
% between the trials
     if length(sequence_time) > 2
         error('Check stimulus times')
     end
     
     %Calculate the Nr of repeats per cycle
     
     Noise_seq_complete = NaN(10*125,Noise_size);
     Repeats_per_epoch = 10;
     
     
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
         
         
     
     
    
     
     
     %% Get only the spikes which are within that time period
        ff = 0;
     f = waitbar(ff,'Stimulus_repeats');  
            
     
     for cc = 1:stx
         ff = cc/stx;
            waitbar(ff);
     
     %Check how many spikes are in the cluster
     temp_spiketimestamps = spiketimestamps(:,cc);
     temp_spiketimestamps(temp_spiketimestamps == 0) = NaN;
     nan_spiketimestamps = isnan(temp_spiketimestamps);
     l_spiketimestamps = length(nan_spiketimestamps) - nnz(nan_spiketimestamps);
     
     %Check which of these spikes accure during the colour noise
     
     colour_begin = Epochs.Epochs_begin(C_epochs);
     colour_end = Epochs.Epochs_end(C_epochs);
     i1 = 1;
     i2 = 0;
     clear true_spikes
     for ee = 1:nr_C_epochs
         
         
         true_spikes_temp = (temp_spiketimestamps >= colour_begin(ee)) .* ...
             (temp_spiketimestamps <= colour_end(ee));
         true_spikes_temp = logical(true_spikes_temp);
         %count how many spikes were found
         i2 = i2 + nnz(true_spikes_temp);
        true_spikes(i1:i2,1) = temp_spiketimestamps(true_spikes_temp);
        true_spikes(i1:i2) = true_spikes(i1:i2) - colour_begin(ee);
        i1 = length(true_spikes) + 1;
        %Set the beginning of the noise stimulus as spiketime 0
        
     end
   
      
   %length of new spiketimestamps matrix 
   tsy = length(true_spikes(:,1));
   
   
   %% Calculate the Kernels
   
   % This variable counts the time since the last frame:
%    
%    F = 0;
%    F1 = 0;
%    stim_idx = 1;
tt = 1;
   for ll = 1:tsy
      ll
      for kk = 1:kernel_bins
          
          spike = true_spikes(ll,1);
          time_frame = spike - (time_window-kk*kbinsize-post_window);
          time_frame = floor(time_frame/time_per_frame);
          
          if time_frame <= 0
              colour_value(kk,(1:Noise_size)) = zeros;
          elseif time_frame > length(Noise_seq_complete(:,1))
               colour_value(kk,(1:Noise_size)) = zeros;
              
              
          else
          colour_value(kk,:) = Noise_seq_complete(time_frame,:);  
          end
      end
      Kernels(:,:,cc) = Kernels(:,:,cc) + colour_value';
      %test_Kernel(:,:,tt) = colour_value;
      tt = tt+1;
      clear colour value
   end
  end
          
              
            
               
              
              
              
              
              
              
              
          end
              
          
          
          
          
          
          
          
       
       
   