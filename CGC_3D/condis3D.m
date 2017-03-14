function d = condis3D(TI,patch,cond)
% The function to calculate the 3D distance

global w_c;

K = ones(size(patch));  K(isnan(patch)) = 0;
patch(isnan(patch)) = 0;
a2 = filter3_v1(K,TI.^2);
ab = filter3_v1(patch,TI).*2;
p = patch.^2;
b2 = sum(p(:));

d1=((a2-ab)+b2);
d2 = max(d1,0);
d = sqrt(d2);
end