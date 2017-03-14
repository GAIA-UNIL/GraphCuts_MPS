function con_dis = dis_con2D_v2(TI,simold,sigforcondi,w)

%% define the filter for unconditioning data
filternon = ones(size(simold));
filternon(isnan(simold)) = 0;
%filternon(isfinite(condi)) = 0;
Nt = sum(filternon(:));

%% define the filter for conditioning data
filterc = isfinite(sigforcondi);
Nc = sum(filterc(:));

%% devide simold to overlap part and conditional part
ovsim = simold;
csim = sigforcondi;
%ovsim(isfinite(condi)) = 0;
ovsim(isnan(ovsim)) = 0;
csim(isnan(sigforcondi)) = 0;

%% the first term of costfunction for unconditioning data
%a2 = filter2(filternon, ovsim.^2, 'valid');
a2 = sum(sum(ovsim.^2));
b2 = filter2(filternon, TI.^2, 'valid');
ab = filter2(ovsim, TI, 'valid')*2;
term1 = a2 + b2 -ab;

if term1 < 0
    disp('negative values occur in the first term of the cost');
end

%% the second term of costfunction for conditioning data
%m2 = filter2(filterc, csim.^2, 'valid');
m2 = sum(sum(csim.^2));
n2 = filter2(filterc, TI.^2, 'valid');
mn = filter2(csim, TI, 'valid')*2;
term2 = m2 + n2 - mn;

if term2 < 0
    disp('negative values occurs in the second term of the cost');
end
%% the patches that exactly match the conditioning data
[pr(:,1),pr(:,2)] = find(term2 == 0);

%% the totalcost for patches in the training image
if Nt == 0 % there is no unconditioning data
    if Nc == 0
        disp('empty region !');
        con_dis = zeros(size(term1));
        % constant value of 0
    else
        con_dis = sqrt(term2)/Nc*w;
    end
else if Nc == 0
        con_dis = sqrt(term1)/Nt*(1-w);
    else
        con_dis = sqrt(term1)/Nt*(1-w) + sqrt(term2)/Nc * w; 
    end
end

end