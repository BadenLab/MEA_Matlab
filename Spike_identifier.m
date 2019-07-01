cluster = 185;
spikes_cluster = cluster_id == cluster;
spike_time = times(spikes_cluster)/Sampling;
spike_shapes = Shapes(spikes_cluster,:);
mean_spikes = mean(Shapes,1);
spike_channel = Channels(spikes_cluster);

%For Channel 13,19
test_channel = Channels == 786;
test_time = times(test_channel);
test_time = test_time / Sampling;
