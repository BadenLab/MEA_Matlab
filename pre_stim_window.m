function Epochs =  pre_stim_window (Epochs, pre_stimulus_window,Mode)


binsize = Epochs.binsize;


additional_bins = ceil(pre_stimulus_window/binsize);
additional_bins = additional_bins*binsize;

if Mode == 0
for ii = 1:Epochs.nr_epochs
    for kk = 1:length(Epochs.stimulus_starts(:,1,1))
Epochs.stimulus_starts_pw(kk,1,ii) = Epochs.stimulus_starts(kk,1,ii) - additional_bins; 
    end
end

elseif Mode == 1
    

for ii = 2:length(Epochs.Epochs_end)
Epochs.stimulus_starts_pw(1,1,ii) = Epochs.Epochs_end(1,(ii-1));
   


end
%find first epochs
    fepoch_idx = Epochs.Epoch_code == single(0.3);
Epochs.stimulus_starts_pw(1,1,fepoch_idx) = Epochs.Epochs_begin(1,fepoch_idx) - additional_bins;

Epochs_begin_pw = squeeze((Epochs.stimulus_starts_pw))';
Epochs.Epochs_duration(1,1,:) = ceil05(Epochs.Epochs_end-Epochs_begin_pw,0.5);
Epochs.Epochs_duration_o(1,1,:) = ceil05(Epochs.Epochs_end-Epochs.Epochs_begin,0.5);



end
