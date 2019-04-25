function out = Find_Epoch (Epochs, Epoch_idx)

[~, ~, Loop_bins_begin, Loop_bins_end] = RT_Loop (Epochs);

%Get the times for the specific epoch

% Loop_time_begin = Loop_time_begin(Epoch_idx);
% Loop_time_end   = Loop_time_end(Epoch_idx);
Loop_bins_begin = Loop_bins_begin(Epoch_idx);
Loop_bins_end   = Loop_bins_end(Epoch_idx);


out = Epochs.Smoothened_averaged_spikes(Loop_bins_begin+1:Loop_bins_end,:);









end


