function [Epochs] = Arduino_Epochs ...
    (Trigger_Ch, Threshold, SamplingFrequency)


diff_trigger = diff(Trigger_Ch);
increase = find(diff_trigger > Threshold);
increase = increase +1;
diff_trigger = double(diff_trigger);

[~,locs] = findpeaks(diff_trigger,'MinPeakHeight',80);

locs = locs-1;
trigger_times_temp = locs/SamplingFrequency;

% trigger_times = (increase-1)/SamplingFrequency;



decrease = find(diff_trigger < -Threshold);
decrease = decrease +1;

stim_time = increase/SamplingFrequency;
stim_timeD = decrease/SamplingFrequency;

stim_time_diff = diff(stim_time);
stim_time_diffD = diff(stim_timeD);

stim_begin = find(stim_time_diff > 0.03);
stim_end = find(stim_time_diffD > 0.03);

% for ii = 1:length(stim_begin)
%    if mod(ii,2) == 0
%        stim_begin(ii) = NaN;
%    else
%        stim_end(ii) = NaN;
%    end
% end
stim_begin = stim_begin +1;


% stim_begin(sbi) = [];
% stim_end(sei) = [];

stim_begin = stim_time(stim_begin);
stim_end = stim_timeD(stim_end);

%Check for the last stimulus end


if stim_end(end) == stim_timeD(end)
    
else
    stim_end(end+1) = stim_timeD(end);
end

start_times(1) = stim_time(1);
sbl = length(stim_begin);
start_times(2:sbl+1) = stim_begin;
end_times = stim_end;


%% Sort the stimuli and assign them to the right stimulus type

cases = length(start_times);

stimulus_begin = NaN(1,(length(start_times)/2));
stimulus_end = NaN(1,(length(end_times)/2));
stimtype_begin = NaN(1,(length(start_times)/2));
stimtype_end = NaN(1,(length(end_times)/2));

for ii = 1:cases
    if mod(ii,2) == 0
        stimulus_begin(ii) = start_times(ii);
        stimulus_end(ii) = end_times(ii);
        stimtype_begin(ii) = NaN;
        stimtype_end(ii) = NaN;
    else
        stimtype_begin(ii) = start_times(ii);
        stimtype_end(ii) = end_times(ii);
        stimulus_begin(ii) = NaN;
        stimulus_end(ii) = NaN;
    end
end

stb = isnan(stimtype_begin);
ste = isnan(stimtype_end);
smb = isnan(stimulus_begin);
sme = isnan(stimulus_end);

stimtype_begin(stb) = [];
stimtype_end(ste) = [];
stimulus_begin(smb) = [];
stimulus_end(sme) = [];

trigger_times = NaN(length(trigger_times_temp(1,:)),length(stimulus_begin));

for ii = 1:length(stimulus_begin)
    trigger_idx1 = (trigger_times_temp > stimulus_begin(ii)).*...
        (trigger_times_temp < stimulus_end(ii));
    trigger_idx1 = logical(trigger_idx1);
    trigger_timesX = trigger_times_temp(trigger_idx1);
    trigger_times(1:(length(trigger_timesX)),ii) = trigger_timesX;
    Nanrow(ii) = nnz(~isnan(trigger_times(:,ii)));
end
Nanrow = max(Nanrow)+1;

trigger_times((Nanrow:end),:) = [];

   
    
    
    

stimcode = ceil05(stimtype_end - stimtype_begin,0.05);

stimtype = NaN(length(stimcode));



trigger_rate = mean(diff(trigger_times));
% trigger_rate = trigger_rate(1,2);
stN = isnan(stimtype);
stimtype(stN) = [];
Epochs.Epochs_begin = stimulus_begin;
Epochs.Epochs_end = stimulus_end;
% Epochs.Epochtext = stimtext;
Epochs.Epoch_code = single(stimcode);
Epochs.trigger_times = trigger_times;
Epochs.trigger_rate = trigger_rate;
Epochs.nr_epochs = length(stimcode);

Epochs.nr_unique_epochs = length(unique(Epochs.Epoch_code));
Epochs.Epochtext = cell(1,Epochs.nr_epochs);
Epochs.nr_stim_repeats = Epochs.nr_epochs/Epochs.nr_unique_epochs;


        
        
    
    

end























