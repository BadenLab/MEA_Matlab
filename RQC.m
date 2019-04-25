function Quality_criteria = RQC (Epochs,Epochs2analyse)
Data = create_trial_responses(Epochs);
Data_size = size(Data);
STx = Data_size(1);
nr_epochs = Data_size(3);

test = exist('Epochs2analyse','var');
if test == 0
    Epochs2analyse = ones(1,nr_epochs);
end


Quality_criteria = zeros(STx,1,nr_epochs);


for ii = 1:STx
    for ss = 1:nr_epochs
        
        Data_temp = squeeze(Data(ii,:,ss,:));
        Data_var = var(Data_temp,0,2,'omitnan');
        mean_of_var = nanmean(Data_var);
        Data_mean = nanmean(Data_temp,2);
        var_of_mean = var(Data_mean,0,1,'omitnan');
         
        
        
        Quality_criteria(ii,1,ss) = var_of_mean/mean_of_var;
    end
end
        

Epochs2analyse = logical(Epochs2analyse);
Quality_criteria = nanmean(Quality_criteria(:,:,Epochs2analyse),3);





end