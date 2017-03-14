function [lineno,flowold,flow,labels,Nadd,check,preseam] = Addseam2D(errori,seamnode,TI, realnew,Startx,Starty,T1,T2)
%ADDSEAM2D scan the seamnode list to find the previous seamnodes covered in
%the new patch. Add the seam nodes to the new overlap. Do cut. 
%   seamnode: the struct list to store the seam nodes:
% seam.dx,seam.dy,seam.dz,seam.cost,seam.adind0,seam.adind1,seam.ind0p0,
% seam.ind1p0,seam.ind0p1,seam.ind1p1

global nbvar
global SGdim

Ngrid = SGdim(1)*SGdim(2);

[ox,oy] = size(errori);
Noverlap = ox*oy;
size1 = size(T1,1);
size2 = size(T2,1);
flowold = 0;

if size(seamnode,1) * size(seamnode,2) == 0
    Nadd = 0;
    lineno = [];
else
    lineno = find(seamnode(:,1)<=Startx+ox-1 & seamnode(:,1)>Startx & seamnode(:,2)<=Starty+oy-1 & seamnode(:,2)>Starty);
    Nadd = size(lineno,1);
end

Esame = edges4connected(ox,oy);
N = Noverlap + Nadd;

Addind0 = seamnode(lineno,5);
Addind1 = seamnode(lineno,6);

preseam = seamnode(lineno,:);
preseam(:,1) = preseam(:,1) - Startx + 1;
preseam(:,2) = preseam(:,2) - Starty + 1;

if Nadd > 0
    [x0,y0] = ind2sub([SGdim(1),SGdim(2)],Addind0);
    [x1,y1] = ind2sub([SGdim(1),SGdim(2)],Addind1);
    tx0=x0-Startx+1; ty0=y0-Starty+1;
    tx1=x1-Startx+1; ty1=y1-Starty+1;
    tind0 = sub2ind([ox,oy],tx0,ty0);
    tind1 = sub2ind([ox,oy],tx1,ty1);
    seamnode_overlap = [seamnode(lineno,1:4),tind0,tind1,seamnode(lineno,7:10)];

    %% Add seam nodes
    Edif = nan(Nadd*4,2);
    Vdif = nan(Nadd*4,1);

    %% connected the seam node to their adjacent nodes
    Edif(1:Nadd,1) = tind0;
    Edif(1:Nadd,2)=[Noverlap+1:N];
    Edif(Nadd+1:2*Nadd,1) = [Noverlap+1:N];
    Edif(Nadd+1:2*Nadd,2) = tind0;
    
    Edif(2*Nadd+1:3*Nadd,1) = [Noverlap+1:N];
    Edif(2*Nadd+1:3*Nadd,2) = tind1;
    Edif(3*Nadd+1:4*Nadd,1) = tind1;
    Edif(3*Nadd+1:4*Nadd,2) = [Noverlap+1:N];
    
    %% assign arc values
    v1 = [];
    v2 = [];
    
    v3 = [];
    v1_new = [];
    v3_new = [];
    for j = 0:nbvar-1
        v1 = [v1,[TI(seamnode_overlap(:,7)+(j*Ngrid));TI(seamnode_overlap(:,8)+(j*Ngrid))]];
        v3 = [v3,[TI(seamnode_overlap(:,9)+(j*Ngrid));TI(seamnode_overlap(:,10)+(j*Ngrid))]];
        v2 = [v2,[realnew(tind0+(j*Noverlap));realnew(tind1+(j*Noverlap))]];
%         v1_new = [v1_new,realnew(tind0+(j*Noverlap))];
%         v3_new = [v3_new,realnew(tind1+(j*Noverlap))];
%         v1 = [v1,TI(seamnode_overlap(:,7)+(j*Ngrid))];
%         v3 = [v3,TI(seamnode_overlap(:,9)+(j*Ngrid))];
    end
    Vdif(1:Nadd) = CutCost(v1,v2);
    Vdif(2*Nadd+1:3*Nadd) = CutCost(v3,v2);
%     Vdif(1:Nadd) = CutCost_v2(v1,v1_new);
%    Vdif(1:Nadd) = abs(seamnode_overlap(:,7)-realnew(tind0)) + abs(seamnode_overlap(:,8)-realnew(tind1));
    Vdif(Nadd+1:2*Nadd) = Vdif(1:Nadd);
%    Vdif(2*Nadd+1:3*Nadd) = abs(seamnode_overlap(:,9)-realnew(tind0))+abs(seamnode_overlap(:,10)-realnew(tind1));
% 	Vdif(2*Nadd+1:3*Nadd) = CutCost_v2(v3,v3_new);
    Vdif(3*Nadd+1:4*Nadd) = Vdif(2*Nadd+1:3*Nadd);
    
    %% remove the old edges
    for uu = 1:Nadd
        removeline0 = find(Esame(:,1)==tind0(uu) & Esame(:,2)==tind1(uu));
        removeline1 = find(Esame(:,1)==tind1(uu) & Esame(:,2)==tind0(uu));
        removeline = [removeline0;removeline1];
        Esame(removeline,:) = [];
    end
    Vsame = errori(Esame(:,1)) + errori(Esame(:,2)) + 9e-9;
    
    E = [Esame;Edif];
    V = [Vsame;Vdif];
    
    T2add = [Noverlap+1:N]';
    ACRV = seamnode_overlap(:,4);
    flowold = sum(ACRV(:));
    T = sparse([T1;T2;T2add],[ones(size1,1);ones(size2,1)*2;ones(Nadd,1)*2],[ones(size1+size2,1)*9e9;ACRV],N,2);
else
    E = Esame;
    V = Vsame;
    T = sparse([T1;T2],[ones(size1,1);ones(size2,1)*2],ones(size1+size2,1)*9e9, N ,2);
end

A = sparse(E(:,1), E(:,2), V, N, N);
[flow,labels] = maxflow(A,T);

if Nadd>0
    check = nan(Nadd,5);
    check(:,1) = labels(Noverlap+1:N);
    check(:,2) = labels(tind0);
    check(:,3) = labels(tind1);
    check(:,4) = tind0;
    check(:,5) = tind1;
    check(:,6) = Vdif(1:Nadd); % cut cost between A and C
    check(:,7) = Vdif(2*Nadd+1:3*Nadd); % cut cost between B and C
end

end
