function Epochs_bined = Epochs_bining_V2(Epochs , spiketimestamps, Mode)
binsize = Epochs.binsize;
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
Epochs_bined = NaN(nr_ST,xST,Nr_bins,Epochs.nr_epochs);
bin_start_1 = NaN(nr_ST,Nr_bins,Epochs.nr_epochs);
bin_end_1 = NaN(nr_ST,Nr_bins,Epochs.nr_epochs);
indexL = NaN(1,xST);
for ii = 1:xST
    indexL(1,ii) = find(spiketimestamps(:,ii),1,'last');
end
    
    
    
    
   







    for bb = 1:Epochs.nr_epochs 
        for ss = 1:nr_ST
%             if isnan(Epochs.stimulus_starts(ss,1,bb)) == 1
%                 break;
%             else
                for kk = 1:Epochs.nr_bins(ss,1,bb)
      
                    bin_start_1(ss,kk,bb) = Epochs.stimulus_starts_pw(ss,1,bb) + binsize * (kk-1);
                    bin_end_1(ss,kk,bb) = bin_start_1(ss,kk,bb) + binsize;
                end
            
        end
    end
    
  
C = ones(1,xST);
B = ones(1,xST);
Bined_epochs_temp = NaN(1,xST);
 gg = 0;
        g = waitbar(gg,'Epochs');

    for bb = 1:Epochs.nr_epochs
       gg = bb/Epochs.nr_epochs;
       waitbar(gg);
        % Go through all the epochs
       ff = 0;
       
%    f = waitbar(ff,'Stimulus_repeats');  

        for ss = 1:nr_ST
%             ff = ss/nr_ST;
%             waitbar(ff);
            %Go through all the stimulus repetition
             
%              if isnan(Epochs.stimulus_starts(ss,1,bb)) == 1
%                 break;
%                 else
                for kk = 1:length(bin_start_1(ss,:,bb))
                    kk
                    if isnan(bin_start_1(ss,kk,bb)) == 1
                        break;
                    else
                        bin_start_temp = bin_start_1(ss,kk,bb);
                        bin_end_temp = bin_end_1(ss,kk,bb);
                   
                    for jj = 1:xST
                        
    %This can be speeded up by only indexing the field with actual values
    
    
    
                            indexa = ((spiketimestamps((C(1,jj):indexL(1,jj)),jj) >= (bin_start_temp))...
                            .* ((spiketimestamps((C(1,jj):indexL(1,jj)),jj) < (bin_end_temp))));
                        %Works
                        
                            Bined_epochs_temp(1,jj) = nnz(indexa); 
                            
                            [~,B(1,jj)] = max(indexa,[],1);
%                             B(1,jj) = B(1,jj) + Bined_epochs_temp(1,jj)-1;
                           
                            B(1,jj) = B(1,jj)-1;
                            C(1,jj) = B(1,jj)+C(1,jj);
                            clear last_index
                            if C(1,jj) >= indexL(1,jj) 
                                C(1,jj) = indexL(1,jj);
                            end
                            end
    
%                             Bined_epochs_temp(1,jj) = nnz((spiketimestamps(:,jj) >= (bin_start_temp))...
%                             .* ((spiketimestamps(:,jj) < (bin_end_temp))));
                            
                            
                            
                            
                            
                            
                    end
                    Epochs_bined(ss,:,kk,bb) = Bined_epochs_temp(1,:);
                    end
                    
                   
%                 end
        end
%              close(f)
        end
    

    close(g)
    
                    

    
       
       



end