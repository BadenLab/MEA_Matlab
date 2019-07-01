


         
time_window = 0.25;
Channel_out = zeros(64,64,round(2*time_window*SamplingFrequency+1));

for a = 1:64
    a
    for b = 1:64
        b
        %Create the name for the channels to be loaded 
        dot = '.';
        Channel_a = string(a);
        Channel_b = string(b);
        text_save = 'Colour_Noise_Ch';
        mat_end = '.mat';
        load_string = strcat(text_save,Channel_a,dot,Channel_b,mat_end);
        Channel_number_1 = sprintf('%02d', a);
        Channel_number_2 = sprintf('%02d', b);
        Channel_str = 'Ch';
        s_s = '_';
        struc_name = 'Data';
        Channel_Field = strcat(Channel_str, Channel_number_1,s_s,Channel_number_2);
        time_window_l = time_window;
      
        
        %exclude Channel 1.1
        test = a*b;
        if test ~= 1
            Data = load(load_string);
            Channel_temp = Data.(Channel_Field);
            length_Channel = length(Channel_temp(1,:));
            SamplingFrequency = Data.SamplingFrequency;
            %Add two time window before and after the beginning of the Channel
            length_Channel_total = round(length_Channel + 2*time_window * SamplingFrequency);
            Channel = zeros(1, length_Channel_total);
            %This adds zeros at the beginning and the end of the channel
            %trace so that we can average also the first and the last
            %spikes
            Channel(1,(time_window_l*round(SamplingFrequency):...
                length_Channel+time_window_l*(round(SamplingFrequency))-1)) = Channel_temp;
            %subtract the begin of the subtimes
            spiketimes_norm = spiketimes - start_time;
            %Check which frame we are looking at
            Frames = spiketimes_norm*SamplingFrequency;
            Frames = round(Frames);
            Frames = Frames + time_window_l*round(SamplingFrequency);
            Sampling_ms = round(SamplingFrequency/1000);
            End_Shift = round(2000*Sampling_ms);
            Channel_out_temp = zeros(length(spiketimes(:,1)),round(2*time_window*SamplingFrequency+1));
             
            parfor s = 1:length(spiketimes(:,1))
                Channel_loop = Channel;
                a1 = a;
                b1 = b;
                
                Begin = Frames(s) - round(time_window_l*SamplingFrequency);
                End = Frames(s) + round(time_window_l*SamplingFrequency);
               
                
                Channel_out_temp1 = Channel_loop(1,Begin:End);
                saturation = Channel_out_temp1(1,:)>500;
                Channel_out_temp1(saturation) = 0;
                
                
                Channel_out_temp(s,:) = Channel_out_temp1;
                                         
            end
            Channel_out(a,b,:) = nansum(Channel_out_temp);
        end
    end
end
         
        
                
                
            
            
%     for s = 1%:length(spiketimes(:,1))
%                 Channel_loop = Channel;
%                 a1 = a;
%                 b1 = b;
%                 s
%                 for f = 1:2500
%                     f1 = -501+f;
%                     Shift = round(f1*Sampling_ms);
%                     Begin = Frames(s) + Shift;
%                     End = Begin + End_Shift;
%                     Channel_out_temp1 = Channel_loop(1,Begin:End);
%                     Channel_out_temp(f,:) = Channel_out_temp(f,:) + Channel_out_temp1;
%                     
%                 end                
%             end         
%             
            
            
        
     
