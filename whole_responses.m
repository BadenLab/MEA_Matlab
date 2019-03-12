function whole_responses = create_trial_responses (Epochs)
%Get the number of cells
STx = gather(length(Epochs.spiketimestamps(1,:)));
Nr_bins = max(Epochs.nr_bins);
%Create a matrix according to how many unique epochs were used
whole_response = NaN(STx,Nr_bins,Epochs.nr_unique_epochs,Epochs.nr_stim_repeat);
Bined_epochs = gather(Epochs.Bined_epochs);

for ii = 1:STx
    ii
    for ee = 1:Epochs.nr_unique_epochs
        %Check which epochs have the same code as the current one
        Epoch_code = Epochs.Epoch_code(ee);
        Epoch_code_w = Epochs.Epoch_code == Epoch_code;
        
    
    
    whole_response_temp = squeeze(Bined_epochs(1,ii,:,...
        (Epoch_code_w)));
    
    whole_response(ii,:,ee,:) = squeeze(whole_response_temp);
        
    
    end
end
    
    


end