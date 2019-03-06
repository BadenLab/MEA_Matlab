function sequence = read_single_stim_seq (Epochs)

for kk = 1:Epochs.nr_unique_epochs
    Epoch_name(kk) = convertCharsToStrings(Epochs.Epochtext(kk));
    indexstr(kk) = 'Epochs.' + Epoch_name(kk);

end
for ii = 1:length(indexstr)
max_seq_temp(ii) = length(eval(indexstr(ii)));
end
max_seq = max(max_seq_temp);
sequence = NaN(max_seq,Epochs.nr_unique_epochs);

for ii = 1:Epochs.nr_unique_epochs
    sequence((1:max_seq_temp(ii)),ii) = eval(indexstr(ii));
end






end