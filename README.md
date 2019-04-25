# MEA_Matlab
This repository contains Matlab scripts for the analysis of sorted spikes (using the HS2 spike sorting) recorded with the MEA


The Version from the 25th of April is modified so that stimuli which contain only one colour noise epoch at the end of the last loop repeat can be analyzed. The change is made to the Responseaverage function. The function tests which epochs are unqique and how often unique epochs are shown (nr_unique). The maximal number of nr_unique defines how often the loop has been repeated. 
At the moment the last epoch gets excluded from the averaging of responses (this shall be changed in future so that it doesnt matter at which position the epoch which is repeated only once is. Today showing the colour noise not at the end of the stimulus will result in an error. 
In addition the function now also returns Epochs.loop_repeats.

The function "create_trial_responses" has been modified in line 6:
whole_response = NaN(STx,Nr_bins,Epochs.nr_unique_epochs,Epochs.loop_repeats);

Also in the plotting section the number of loop repeats is no taken from that variable
stimulus_repeat = Epochs.loop_repeats;




