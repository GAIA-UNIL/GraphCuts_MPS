function [filtermap, Noverlap] = distance2D_v1(TI,simold)
% the distance between selected patch and the training image
Noverlap = size(simold,1) * size(simold,2);
filter = ones(size(simold));
filter(isnan(simold)) = 0;
fsim = simold; fsim(isnan(simold)) = 0;

a2 = filter2(filter,fsim.^2,'valid');
b2 = filter2(filter,TI.^2,'valid');
ab = filter2(fsim,TI,'valid')*2;

term = a2+b2-ab;
if term < 0
    disp('negative values occur in the first term of the cost');
end
filtermap = sqrt(term)/Noverlap;
end

