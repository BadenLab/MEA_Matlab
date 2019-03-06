function Epochs = Bin_fitting_V2 (Epochs, Mode)

% for ii = 1:Epochs.nr_epochs
%     for kk = 1:length(Epochs.stimulus_starts_pw(:,1,1))
% Epochs.nr_bins(kk,1,ii) = ceil((Epochs.stimulus_ends(kk,:,ii) - Epochs.stimulus_starts(kk,:,ii))...
%     /binsize);
%     end
binsize = Epochs.binsize;
if Mode == 0
    for ii = 1:Epochs.nr_epochs
    for kk = 1:length(Epochs.stimulus_starts_pw(:,1,1))
Epochs.nr_bins(kk,1,ii) = ceil((Epochs.stimulus_ends_pw(kk,:,ii) - Epochs.stimulus_starts_pw(kk,:,ii))...
    /binsize);
    end
    end
elseif Mode == 1
    
    for ii = 1:Epochs.nr_epochs
%         Epochs.nr_bins(1,1,ii) = ceil((Epochs.Epochs_end(1,ii) - Epochs.stimulus_starts_pw(1,1,ii))...
%     /binsize);
% 
 Epochs.nr_bins(1,1,ii) = Epochs.Epochs_duration(1,1,ii)/binsize;
    end
end
    
   

end