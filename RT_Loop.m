function [Loop_time_begin, Loop_time_end, Loop_bins_begin, Loop_bins_end] = RT_Loop (Epochs)

% This function calculates the beginning and end times and bins for the
% averaged whole loop


Loop_time_begin = zeros(1,Epochs.nr_epochs_plot);
Loop_time_end = zeros(1,Epochs.nr_epochs_plot);

Epochs_begin_temp = Epochs.stimulus_starts_pw(1,1:Epochs.nr_epochs_plot);
Epochs_end_temp = Epochs.Epochs_end(1,1:Epochs.nr_epochs_plot);

Loop_time_begin = round(Epochs_begin_temp - Epochs_begin_temp(1));
Loop_time_end = round(Epochs_end_temp - Epochs_begin_temp(1));



Loop_bins_begin = Loop_time_begin./Epochs.binsize;
Loop_bins_end = Loop_time_end./Epochs.binsize;




end