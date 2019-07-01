%% Quality control of Kernels

function [KQI,KQI_C,KQI_past] = Quality_Kernels(Kernels,colours,border)

Kernels = gather(Kernels);
%how many cells?
l_K = length(Kernels(1,1,:));
%how many bins
% b_K = length(Kernels(1,:,1));



KQI_temp = zeros(l_K,colours);
KQI_size = zeros(l_K,colours);
%loop over cells
for kk = 1:l_K
    

%loop over different colours
for ii = 1:colours
    
    Kernels_temp = Kernels(ii,:,kk); 
    
    %find the minmal and maximal stim count
    min_temp = min(Kernels_temp);
    max_temp = max(Kernels_temp);
    
    %check if min is zeros
    if min_temp == 0
        difference_stim = 0;
    else
    
    difference_stim = max_temp/min_temp - 1;
    end
    
    %if the difference is smaller than 10% the Kernels will be counted as
    %not good
    if difference_stim < border
        continue
        
    else
        KQI_size(kk,ii) = difference_stim;
        
        %If the difference is bigger than 10% we apply the next test
        KQI_temp(kk,ii) = std(Kernels_temp);
        
        
        
        
        
    end
    
    
end
        
        
    
    
    




end

KQI_C = KQI_size .* KQI_temp;
bar(KQI_C,'DisplayName','KQI_C')

KQI = max(KQI_C,[],2);
KQI_past = KQI > (nanmean(KQI)+std(KQI));
KQI_past = find(KQI_past)';



end