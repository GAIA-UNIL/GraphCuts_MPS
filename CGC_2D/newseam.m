function [newseamnode] = newseam(AC,startx,starty,startz,labelnew, labelold, errori,TI)
%NEWSEAM the new seam nodes identified from the cut label
%   newseamnode
% x,y,z,cost,ind0,ind1,TI_ind0,TI_ind0_adjacent, TI_ind1, TI_ind1_adjacent

[dimx,dimy,dimz] = size(errori);
global SGdim
global nbvar
Ngrid = SGdim(1)*SGdim(2)*SGdim(3);

ind0 = AC(:,2); ind1 = AC(:,3);
[x0,y0,z0] = ind2sub([dimx,dimy,dimz],ind0);
[x1,y1,z1] = ind2sub([dimx,dimy,dimz],ind1);

tx0 = x0 + startx - 1;
ty0 = y0 + starty - 1;
tz0 = z0 + startz - 1;

tx1 = x1 + startx - 1;
ty1 = y1 + starty - 1;
tz1 = z1 + startz - 1;

tx = 0.5*(tx0+tx1);
ty = 0.5*(ty0+ty1);
tz = 0.5*(tz0+tz1);

%% 
cutcost = errori(ind0) + errori(ind1);

%% TI

tind0 = sub2ind(SGdim(1:3),tx0,ty0,tz0);
tind1 = sub2ind(SGdim(1:3),tx1,ty1,tz1);

v1=[]; v2=[];
for k = 0:nbvar-1
    v1 = [v1,[TI(labelold(ind0)+k*Ngrid);TI(labelold(ind1)+k*Ngrid)]];
    v2 = [v2,[TI(labelnew(ind0)+k*Ngrid);TI(labelnew(ind1)+k*Ngrid)]];
end
cutcost = CutCost(v1,v2);

newseamnode = [tx,ty,tz,cutcost,tind0,tind1,labelold(ind0),labelold(ind1),labelnew(ind0),labelnew(ind1)];
% notice that the columns 9 and 10 are connected but the columns 7 and 8
% maybe not connected
end

