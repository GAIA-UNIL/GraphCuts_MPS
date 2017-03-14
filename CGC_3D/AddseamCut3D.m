function [oldseamnode,newseamnode,label,preseam] = AddseamCut3D(errori,seamnode,TI, realnew,labelold, labelnew,Startx,Starty,Startz,T1,T2)

% the function to make a new cut taking account of previous seam nodes
% seamnode: the list of seam nodes [x,y,z,cutcost,]

%global nbvar
global SGdim
global TIdim
Ngrid = TIdim(1)*TIdim(2)*TIdim(3);

% global variabletype
% global w_v

[ox,oy,oz] = size(errori);
Noverlap = ox*oy*oz;
size1 = size(T1,1);
size2 = size(T2,1);

oldseamnode = seamnode;

if size(seamnode,1) * size(seamnode,2) == 0
    Nadd = 0;
    lineno = [];
    preseam = [];
else
    lineno = find(seamnode(:,1)<=Startx+ox-1 & seamnode(:,1)>=Startx & seamnode(:,2)<=Starty+oy-1 & seamnode(:,2)>=Starty & seamnode(:,3)<=Startz+oz-1 & seamnode(:,3)>=Startz);
    Nadd = size(lineno,1)*size(lineno,2);
    Addind0 = seamnode(lineno,5); 
    Addind1 = seamnode(lineno,6);
end

%Esame = edges4connected(ox,oy);
Esame = edges6connected3D(realnew);

if Nadd > 0
    [x0,y0,z0] = ind2sub([SGdim(1),SGdim(2),SGdim(3)],Addind0);
    [x1,y1,z1] = ind2sub([SGdim(1),SGdim(2),SGdim(3)],Addind1);
    tx0=x0-Startx+1; ty0=y0-Starty+1; tz0=z0-Startz+1;
    tx1=x1-Startx+1; ty1=y1-Starty+1; tz1=z1-Startz+1;
    tind0 = sub2ind([ox,oy,oz],tx0,ty0,tz0);
    tind1 = sub2ind([ox,oy,oz],tx1,ty1,tz1);

    N = Noverlap + Nadd;
    seamnode_overlap = seamnode(lineno,:);
    seamnode(lineno,:) = [];
%    seamnode_overlap = [seamnode(lineno,1:4),tind0,tind1,seamnode(lineno,7:10)];
    pre_x = 0.5*(tx0+tx1);
    pre_y = 0.5*(ty0+ty1);
    pre_z = 0.5*(tz0+tz1);
    preseam = [pre_x,pre_y,pre_z,seamnode_overlap(:,4)];
    %% Add seam nodes
    Edif = nan(Nadd*4,2);
    Vdif = nan(Nadd*4,1);
    
    %% connected the seam node to their adjacent nodes
    Edif(1:Nadd,1) = tind0;
    Edif(1:Nadd,2)=[Noverlap+1:N]';
    Edif(Nadd+1:2*Nadd,1) = [Noverlap+1:N]';
    Edif(Nadd+1:2*Nadd,2) = tind0;
    
    Edif(2*Nadd+1:3*Nadd,1) = [Noverlap+1:N]';
    Edif(2*Nadd+1:3*Nadd,2) = tind1;
    Edif(3*Nadd+1:4*Nadd,1) = tind1;
    Edif(3*Nadd+1:4*Nadd,2) = [Noverlap+1:N]';
    
    %% assign arc values
    if sum(sum(round(seamnode_overlap(:,7:10))==seamnode_overlap(:,7:10))) ~= 4*Nadd
        disp('check seamnode index');
    end
    if sum(sum(seamnode_overlap(:,7:10)<1)) > 0
        disp('check seamnode index');
    end
       
    %% check wether realnew(tind0+(j*Noverlap)) == TI
    costAC = abs(TI(seamnode_overlap(:,7))-realnew(tind0))+abs(TI(seamnode_overlap(:,8))-realnew(tind1))+eps;
    costBC = abs(TI(seamnode_overlap(:,9))-realnew(tind0))+abs(TI(seamnode_overlap(:,10))-realnew(tind1))+eps;
    costAB = seamnode_overlap(:,4);

    Vdif(1:Nadd) = costAC;
    Vdif(Nadd+1:2*Nadd) = Vdif(1:Nadd);
    
    Vdif(2*Nadd+1:3*Nadd) = costBC;
    Vdif(3*Nadd+1:4*Nadd) = Vdif(2*Nadd+1:3*Nadd);

    %% remove the old edges
    for uu = 1:Nadd
        removeline0 = find(Esame(:,1)==tind0(uu) & Esame(:,2)==tind1(uu));
        removeline1 = find(Esame(:,1)==tind1(uu) & Esame(:,2)==tind0(uu));
        removeline = [removeline0;removeline1];
        Esame(removeline,:) = [];
    end
    
    Vsame = errori(Esame(:,1)) + errori(Esame(:,2))+ eps;
    
    E = [Esame;Edif];
    V = [Vsame;Vdif];
    
    T2add = [Noverlap+1:N]';
    ACRV = seamnode_overlap(:,4);
%     if sum(ACRV == costAB) ~= Nadd
%         disp('check cost'); % the difference is 4*eps for each unequal element
%     end
    flowold = sum(ACRV(:));
    T = sparse([T1;T2;T2add],[ones(size1,1);ones(size2,1)*2;ones(Nadd,1)*2],[ones(size1+size2,1)*9e9;ACRV],N,2);
else
    N = Noverlap;
    E = Esame;
    V = errori(E(:,1)) + errori(E(:,2))+eps;
    T = sparse([T1;T2],[ones(size1,1);ones(size2,1)*2],ones(size1+size2,1)*9e9, N ,2);
    preseam = [];
end

A = sparse(E(:,1), E(:,2), V, N, N);
[flow,labels] = maxflow(A,T);

%% check the new cut
label = reshape(labels(1:Noverlap),ox,oy,oz);
%[AC,Ne] = Addconnected2D(label);
AC = Addconnected3D(label);

% figure(21); ViewGrid(label);
% axis equal tight; view(30,30)

%% get the new seam nodes
new_ind0 = AC(:,2); new_ind1 = AC(:,3);
[new_x0,new_y0,new_z0] = ind2sub([ox,oy,oz],new_ind0);
[new_x1,new_y1,new_z1] = ind2sub([ox,oy,oz],new_ind1);

t_new_x0 = new_x0 + Startx - 1;
t_new_y0 = new_y0 + Starty - 1;
t_new_z0 = new_z0 + Startz - 1;

t_new_x1 = new_x1 + Startx - 1;
t_new_y1 = new_y1 + Starty - 1;
t_new_z1 = new_z1 + Startz - 1;

new_sx = (t_new_x0 + t_new_x1)*0.5;
new_sy = (t_new_y0 + t_new_y1)*0.5;
new_sz = (t_new_z0 + t_new_z1)*0.5;

t_new_ind0 = sub2ind(SGdim(1:3),t_new_x0, t_new_y0, t_new_z0);
t_new_ind1 = sub2ind(SGdim(1:3),t_new_x1, t_new_y1, t_new_z1);
cutcost = errori(new_ind0) + errori(new_ind1) + eps;

%labelold(isnan(labelold)) = labelnew(isnan(labelold));
% A_ind =[labelold(new_ind0),labelold(new_ind1)];
% [nux,nuy] = find(isnan(A_ind));
% if size(nux,1) ~= 0
%     w_new_x0 = new_x0(nux);
%     w_new_y0 = new_y0(nux);
%     w_new_x1 = new_x1(nux);
%     w_new_y1 = new_y1(nux);
% end

oldpatchindex0 = labelold(new_ind0); oldpatchindex1 = labelold(new_ind1);
newpatchindex0 = labelnew(new_ind0); newpatchindex1 = labelnew(new_ind1);

if sum(sum(round(oldpatchindex0)==oldpatchindex0)) ~= size(new_ind0)
    disp('check new seam node index');
end
if sum(sum(round(oldpatchindex1)==oldpatchindex1)) ~= size(new_ind1)
    disp('check new seam node index');
end
if sum(sum(round(newpatchindex0)==newpatchindex0)) ~= size(new_ind0)
    disp('check new seam node index');
end
if sum(sum(round(newpatchindex1)==newpatchindex1)) ~= size(new_ind1)
    disp('check new seam node index');
end

seam_new = [new_sx,new_sy,new_sz,cutcost,t_new_ind0,t_new_ind1,labelold(new_ind0),labelold(new_ind1),labelnew(new_ind0),labelnew(new_ind1)];


%keep_old = [];
keep_old1 = [];

if Nadd > 0
    % labels of the seam node, ind0 and ind1
    check = [labels(Noverlap+1:N),labels(tind0),labels(tind1)];
    
    %% a new way to check the old seam
    
    for nn = 1:Nadd
        if check(nn,:) == [0 0 0];
            % keep the old cut;
            keep_old1 = [keep_old1;seamnode_overlap(nn,:)];
        else if check(nn,:) == [1,1,0]
                % combination of BC
                line = find(AC(:,3)==tind0(nn) & AC(:,2)==tind1(nn));
%                 newseamnodeindeforoldpatch = [oldpatchindex0(line),oldpatchindex1(line),newpatchindex0(line),newpatchindex1(line)];
%                 oldseamnodeindeforoldpatch = seamnode_overlap(nn,7:10);
                % zhe liang ge zhi ying gai shi fan de
                
                seam_new(line,4) = costBC(nn);
                seam_new(line,8) = seamnode_overlap(nn,9);
            else if check(nn,:) == [1 0 1]
                     line = find(AC(:,2)==tind0(nn) & AC(:,3)==tind1(nn));
%                     newseamnodeindeforoldpatch = [oldpatchindex0(line),oldpatchindex1(line),newpatchindex0(line),newpatchindex1(line)];
%                     oldseamnodeindeforoldpatch = seamnode_overlap(nn,7:10);
                    
                    seam_new(line,4) = costAC(nn);
                    seam_new(line,8) = seamnode_overlap(nn,8);
                else if check(nn,:) ~= [1 1 1]
                        disp('unexpect results');
                        pause;
                    end
                end
            end
        end
    end % loop for each seam node 


end

%newseamnode = seamnode;
newseamnode = [seamnode;seam_new;keep_old1];

costmap = zeros(SGdim(1:3));
costmap(newseamnode(:,5)) = 0.5 * newseamnode(:,4);
costmap(newseamnode(:,6)) = 0.5 * newseamnode(:,4);

% figure(22); clf;
% ViewGrid(costmap); axis equal tight; view(30,30);

end