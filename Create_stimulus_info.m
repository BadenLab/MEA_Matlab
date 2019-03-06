function [Epochs] = Create_stimulus_info (Epochs,fff_repeat, r_repeat)

stimulus_matrix = max(fff_repeat,r_repeat);
stimulus_starts = NaN(stimulus_matrix,1,length(Epochs.nr_epochs));
stimulus_ends = NaN(stimulus_matrix,1,length(Epochs.nr_epochs));


for aa = 1:length(Epochs.Epoch_code)
    
    
    
    switch Epochs.Epoch_code(aa)
        case {0.5, 1, 1.5, 2, 2.5, 3}
           
            Epochs.stimulus_time(aa) = (Epochs.Epochs_end(aa)-Epochs.Epochs_begin(aa))...
                /fff_repeat;
           
            for ii = 1:fff_repeat
                stimulus_starts(ii,1,aa) = Epochs.Epochs_begin(aa)...
                    +(ii-1)*Epochs.stimulus_time(aa);
                stimulus_ends(ii,1,aa) = stimulus_starts(ii,1,aa) + Epochs.stimulus_time(aa);
            end
           
        case 3.5
            Epochs.stimulus_time(aa) = (Epochs.Epochs_end(aa)-Epochs.Epochs_begin(aa))...
                /r_repeat;
            
           
            
            for ii = 1:r_repeat
                stimulus_starts(ii,1,aa) = Epochs.Epochs_begin(aa)...
                    +(ii-1)*Epochs.stimulus_time(aa);
                stimulus_ends(ii,1,aa) = stimulus_starts(ii,1,aa) + Epochs.stimulus_time(aa);
            end
            
        case 4
            stimulus_starts(1,1,aa) = Epochs.Epochs_begin(aa);
            stimulus_ends(1,1,aa) = Epochs.Epochs_end(aa);
            
           
            
            
            
            
            
        
    end
   
    
    
    
    
end
%This is commented out because it will be calculated later
% stimulus_starts_pw = stimulus_starts(:,:,:) - 1;
stimulus_ends_pw = stimulus_ends(:,:,:) - 1;

Epochs.stimulus_starts = stimulus_starts;
% Epochs.stimulus_starts_pw = stimulus_starts_pw;
Epochs.stimulus_ends = stimulus_ends;
Epochs.stimulus_ends_pw = stimulus_ends_pw;

end