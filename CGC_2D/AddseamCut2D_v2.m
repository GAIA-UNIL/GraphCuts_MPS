function [oldseamnode,newseamnode,label,preseam] = AddseamCut2D_v2(errori,seamnode,TI, realnew,labelold, labelnew,Startx,Starty,T1,T2)

% the function to make a new cut taking account of previous seam nodes
% seamnode: the list of seam nodes [x,y,z,cutcost,]

global nbvar
global TIdim
global SGdim
Ngrid = TIdim(1)*TIdim(2);
[ox,oy] = size(errori);
Noverlap = ox*oy;
size1 = size(T1,1);
size2 = size(T2,1);

oldseamnode = seamnode;

if size(seamnode,1) * size(seamnode,2) == 0
    Nadd = 0;
    lineno = [];
    preseam = [];
else
    lineno = find(seamnode(:,1)<=Startx+ox-1 & seamnode(:,1)>=Startx & seamnode(:,2)<=Starty+oy-1 & seamnode(:,2)>=Starty);
    Nadd = size(lineno,1);
    Addind0 = seamnode(lineno,5);
    Addind1 = seamnode(lineno,6);
end

Esame = edges4connected(ox,oy);
N = Noverlap + Nadd;



if Nadd > 0
    [x0,y0] = ind2sub([SGdim(1),SGdim(2)],Addind0);
    [x1,y1] = ind2sub([SGdim(1),SGdim(2)],Addind1);
    tx0=x0-Startx+1; ty0=y0-Starty+1;
    tx1=x1-Startx+1; ty1=y1-Starty+1;
    tind0 = sub2ind([ox,oy],tx0,ty0);
    tind1 = sub2ind([ox,oy],tx1,ty1);
    seamnode_overlap = seamnode(lineno,:);
%    seamnode_overlap = [seamnode(lineno,1:4),tind0,tind1,seamnode(lineno,7:10)];
    pre_x = 0.5*(tx0+tx1);
    pre_y = 0.5*(ty0+ty1);
    pre_z = ones(Nadd,1);
    preseam = [pre_x,pre_y,pre_z,seamnode_overlap(:,4)];
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
    
    Edif(4*Nadd+1:5*Nadd,1) = [Noverlap+1:N];
    Edif(4*Nadd+1:5*Nadd,2) = T2(1);
    Edif(5*Nadd+1:6*Nadd,2) = [Noverlap+1:N];
    Edif(5*Nadd+1:6*Nadd,1) = T2(1);
    %% assign arc values
    v1 = [];
    v2 = [];
    v3 = [];
    
    for j = 0:nbvar-1
        v1 = [v1,[TI(seamnode_overlap(:,7)+(j*Ngrid));TI(seamnode_overlap(:,8)+(j*Ngrid))]];
        v3 = [v3,[TI(seamnode_overlap(:,9)+(j*Ngrid));TI(seamnode_overlap(:,10)+(j*Ngrid))]];
        v2 = [v2,[realnew(tind0+(j*Noverlap));realnew(tind1+(j*Noverlap))]];
    end
    costAC = CutCost(v1,v2);
    costBC = CutCost(v3,v2);
    costAB = CutCost(v1,v3);
    Vdif(1:Nadd) = costAC;
    Vdif(2*Nadd+1:3*Nadd) = costBC;
    Vdif(Nadd+1:2*Nadd) = Vdif(1:Nadd);
    Vdif(3*Nadd+1:4*Nadd) = Vdif(2*Nadd+1:3*Nadd);
    ACRV = seamnode_overlap(:,4);
    Vdif(4*Nadd+1:6*Nadd) = [ACRV;ACRV];

    
    %% remove the old edges
    for uu = 1:Nadd
        removeline0 = find(Esame(:,1)==tind0(uu) & Esame(:,2)==tind1(uu));
        removeline1 = find(Esame(:,1)==tind1(uu) & Esame(:,2)==tind0(uu));
        removeline = [removeline0;removeline1];
        Esame(removeline,:) = [];
    end
    
    Vsame = errori(Esame(:,1)) + errori(Esame(:,2))+9e-9;
    
    E = [Esame;Edif];
    V = [Vsame;Vdif];
    
%     T2add = [Noverlap+1:N]';
%     ACRV = seamnode_overlap(:,4);
    flowold = sum(ACRV(:));
%    T = sparse([T1;T2;T2add],[ones(size1,1);ones(size2,1)*2;ones(Nadd,1)*2],[ones(size1+size2,1)*9e9;ACRV],N,2);
else
    E = Esame;
    V = errori(E(:,1)) + errori(E(:,2))+9e-9;
    preseam = [];
%    T = sparse([T1;T2],[ones(size1,1);ones(size2,1)*2],ones(size1+size2,1)*9e9, N ,2);
end

A = sparse(E(:,1), E(:,2), V, N, N);
T = sparse([T1;T2],[ones(size1,1);ones(size2,1)*2],ones(size1+size2,1)*9e9, N ,2);
[flow,labels] = maxflow(A,T);

%% check the new cut
label = reshape(labels(1:Noverlap),ox,oy);
[AC,Ne] = Addconnected2D(label);

%% get the new seam nodes
new_ind0 = AC(:,2); new_ind1 = AC(:,3);
[new_x0,new_y0] = ind2sub([ox,oy],new_ind0);
[new_x1,new_y1] = ind2sub([ox,oy],new_ind1);

t_new_x0 = new_x0 + Startx - 1;
t_new_y0 = new_y0 + Starty - 1;

t_new_x1 = new_x1 + Startx - 1;
t_new_y1 = new_y1 + Starty - 1;

new_sx = (t_new_x0 + t_new_x1)*0.5;
new_sy = (t_new_y0 + t_new_y1)*0.5;
new_sz = ones(size(new_sx,1),1);

t_new_ind0 = sub2ind(SGdim(1:2),t_new_x0, t_new_y0);
t_new_ind1 = sub2ind(SGdim(1:2),t_new_x1, t_new_y1);
cutcost = errori(new_ind0) + errori(new_ind1) + 9e-9;

seam_new = [new_sx,new_sy,new_sz,cutcost,t_new_ind0,t_new_ind1,labelold(new_ind0),labelold(new_ind1),labelnew(new_ind0),labelnew(new_ind1)];

keep_old = [];

if Nadd > 0
    % labels of the seam node, ind0 and ind1
    check = [labels(Noverlap+1:N),labels(tind0),labels(tind1)];
%     costAC = Vdif(1:Nadd);
%     costBC = Vdif(2*Nadd+1:3*Nadd);
    
    for nn = 1:Nadd
        ch = [num2str(check(nn,1)),num2str(check(nn,2)),num2str(check(nn,3))];
        switch ch
            case '000'
                % keep the old cut
                keep_old = [keep_old;seamnode_overlap(nn,:)];
            case '110'
                % make a new cut at the same place
                % remove the old one and keep the new one
 %               new_i = seamnode_overlap(nn,:);
                changeline = find(seam_new(:,5)==seamnode_overlap(nn,6) & seam_new(:,6)==seamnode_overlap(nn,5));
                seam_new(changeline,4) = costBC(nn);
                seam_new(changeline,8) = seamnode_overlap(nn,9);
            case '101'
                % make a new cut at the same place
                % remove the old one and keep the new one
                changeline = find(seam_new(:,5)==seamnode_overlap(nn,5) & seam_new(:,6)==seamnode_overlap(nn,6));
                seam_new(changeline,4) = costAC(nn);
                seam_new(changeline,8) = seamnode_overlap(nn,8);
            case '111'
                % remove the old cut
            otherwise
                disp('unexpected cut!'); 
                disp(ch);
        end
    end

end

newseamnode = seamnode;
newseamnode(lineno,:) = [];
newseamnode = [newseamnode;seam_new;keep_old];

end