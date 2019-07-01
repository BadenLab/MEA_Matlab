function Epochs_bined = Epochs_bining_V4(m, Epochs, spiketimestamps, Mode)

%%
% yST = length(spiketimestamps(:,1));
xST = length(spiketimestamps(1,:));
if Mode == 0
    
    nr_ST = length(Epochs.stimulus_starts(:,1,1));
elseif Mode == 1
    nr_ST = 1;
%     Epochs.stimulus_starts_pw = NaN(1,1,length(Epochs.Epochs_begin(1,:)));
%     Epochs.stimulus_starts_pw(1,1,:) = Epochs.Epochs_begin;
end
Nr_bins = max(max(Epochs.nr_bins,[],3));
m.Epochs_bined(nr_ST,xST,Nr_bins,Epochs.nr_epochs) = NaN;

indexL = NaN(1,xST);
for ii = 1:xST
    indexL(1,ii) = find(spiketimestamps(:,ii),1,'last');
end

bin_start_1 = NaN(1,Nr_bins,Epochs.nr_epochs);
bin_end_1 = NaN(1,Nr_bins,Epochs.nr_epochs);

%First we need to calculate the edges for the histocount function 
for bb = 1:Epochs.nr_epochs 
        for ss = 1:nr_ST
%             if isnan(Epochs.stimulus_starts(ss,1,bb)) == 1
%                 break;
%             else
                for kk = 1:Epochs.nr_bins(ss,1,bb)
      
                    bin_start_1(ss,kk,bb) = Epochs.stimulus_starts_pw(ss,1,bb) + Epochs.binsize * (kk-1);
                    bin_end_1(ss,kk,bb) = bin_start_1(ss,kk,bb) + Epochs.binsize;
                end
            
        end
end
%This is stupid but in the histocount function the first edge has to be the
%left edge of the first bin but the last edge has to be the right edge of
%the last bin. So we have to change the last entry in the bin_start_1
%matrix accordingly

for bb = 1:Epochs.nr_epochs
bin_start_1(1,Epochs.nr_bins(1,1,bb),bb) = bin_end_1(1,Epochs.nr_bins(1,1,bb),bb);
end



    
 gg = 0;
        g = waitbar(gg,'Epochs');

for bb = 1:Epochs.nr_epochs
     gg = bb/Epochs.nr_epochs;
       waitbar(gg);
     for kk = 1:xST
        %Take one cell out of the spiketimestamps
        spikes_temp = full(spiketimestamps((1:indexL(1,kk)),kk));
        %Check which spikes fall into the time of the stimulus
        spikes_temp_time = (spikes_temp>Epochs.stimulus_starts_pw(1,1,bb)).*...
            (spikes_temp<Epochs.Epochs_end(1,bb));
        spikes_temp_time = logical(spikes_temp_time);
        %Consider only spikes which fall into the time of the stimulus
        spikes_temp = spikes_temp(spikes_temp_time);
        %We have to fill up the matrix as the Epochs_bined matrix is
        %1:xST:(max(nr_bins):...
        [Bined_temp edges] = histcounts(spikes_temp,bin_start_1(1,(1:Epochs.nr_bins(1,1,bb)),bb));
        last_Bined_temp = length(Bined_temp(1,:));
        Bined_temp(1,(last_Bined_temp+1:Nr_bins)) = NaN;
        m.Epochs_bined(1,kk,:,bb) = Bined_temp;
        test1 = nnz(spikes_temp)
        test2 = sum(Bined_temp)
        spikes_size = length(spikes_temp(:,1));
        test_traces(bb,1:spikes_size) = spikes_temp;

       
    end
end
       close(g) 
        
        
    
   



    
                    

    
       
       



end