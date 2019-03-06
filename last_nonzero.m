function [last_nonz] = last_nonzero(test_for_zero,nr_repeats)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
last_nonz = zeros;
for ii = 1:nr_repeats
[r,c] = find(test_for_zero(:,ii),1,'last');
if last_nonz > r 
    continue
else
    last_nonz = r;
    
end
end



end 

