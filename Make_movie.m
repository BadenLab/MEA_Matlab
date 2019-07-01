%smoothing raw data
Channel_out_s = zeros(64,64,length(Channel_out(1,1,:)));
parfor a = 1:64
    for b = 1:64
        Data = squeeze(Channel_out(a,b,:));
        Data = smooth(Data,1000);
        if mean(Data) > 10000
            Data = zeros(length(Data(:,1)),1);
        end
        Data(1:200,1) = mean(Data);
        Data(end-200:end,1) = mean(Data);
        Channel_out_s(a,b,:) = Data
        
    end
end

Channel_out_n = normalize(Channel_out_s,3);

I = mat2gray(Channel_out_n);
% Channel_out_ex = Channel_out;
% Channel_out_ex(Channel_out_ex>10000) = 0;
% Channel_out_s = Channel_
%normalize
%make movie gpu array

%Change the value of channels which are saturated the whole time to 0.5
for a = 1:64
    for b = 1:64
        for_test = I(a,b,:);
        if all(for_test == 1)
            I(a,b,:) = 0.5;
        end
    end
end

%Show the movie
handle = implay(I,100);
%Define colourmap
map = hot(255);
%Load colour map into the movie
%max = black
handle.Visual.ColorMap.Map(1:255,:) = map;




