function [oldseamnode,newseamnode,label,preseam] = AddseamCut2D(errori,seamnode,TI, realnew,labelold, labelnew,Startx,Starty,T1,T2)

% the function to make a new cut taking account of previous seam nodes
% seamnode: the list of seam nodes [x,y,z,cutcost,]

global nbvar
global SGdim
global TIdim
Ngrid = TIdim(1)*TIdim(2);

% global variabletype
% global w_v

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

if Nadd > 0
    [x0,y0] = ind2sub([SGdim(1),SGdim(2)],Addind0);
    [x1,y1] = ind2sub([SGdim(1),SGdim(2)],Addind1);
    tx0=x0-Startx+1; ty0=y0-Starty+1;
    tx1=x1-Startx+1; ty1=y1-Starty+1;
    tind0 = sub2ind([ox,oy],tx0,ty0);
    tind1 = sub2ind([ox,oy],tx1,ty1);
    
    %% remove the line on the boundary
%     boun = [T1;T2];
%     boun_tind0 = intersect(tind0,boun);
%     boun_tind1 = intersect(tind1,boun);
%     n0=size(boun_tind0,1)*size(boun_tind0,2);
%     n1=size(boun_tind1,1)*size(boun_tind1,2);
%     reline0 = []; reline1 = [];
%     if n0 ~=0
%         for k = 1:n0
%             re = find(tind0 == boun_tind0(k));
%             reline0 = [reline0;re];
%         end
%     end
%     if n1 ~= 0
%         for k = 1:n1
%             re = find(tind1 == boun_tind1(k));
%             reline1 = [reline1;re];
%         end
%     end
%     reline = union(reline0,reline1);
%     sizeremove = size(reline,1);
%     if sizeremove > 0
%         lineno(reline) = [];
%         Addind0(reline) = []; Addind1(reline) = [];
%         tx0(reline)=[]; tx1(reline)=[];
%         ty0(reline)=[]; ty1(reline)=[];
%         tind0(reline)=[]; tind1(reline)=[];
%         Nadd = Nadd - sizeremove;
%     end
end

if Nadd > 0    
    N = Noverlap + Nadd;
    seamnode_overlap = seamnode(lineno,:);
    seamnode(lineno,:) = [];
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
    Edif(1:Nadd,2)=[Noverlap+1:N]';
    Edif(Nadd+1:2*Nadd,1) = [Noverlap+1:N]';
    Edif(Nadd+1:2*Nadd,2) = tind0;
    
    Edif(2*Nadd+1:3*Nadd,1) = [Noverlap+1:N]';
    Edif(2*Nadd+1:3*Nadd,2) = tind1;
    Edif(3*Nadd+1:4*Nadd,1) = tind1;
    Edif(3*Nadd+1:4*Nadd,2) = [Noverlap+1:N]';
    
    %% assign arc values
    v1 = [];
    v2 = [];
    v3 = [];
    
    if sum(sum(round(seamnode_overlap(:,7:10))==seamnode_overlap(:,7:10))) ~= 4*Nadd
        disp('check seamnode index');
    end
    if sum(sum(seamnode_overlap(:,7:10)<1)) > 0
        disp('check seamnode index');
    end
    
    for j = 0:nbvar-1
        v1 = [v1,[TI(seamnode_overlap(:,7)+(j*Ngrid));TI(seamnode_overlap(:,8)+(j*Ngrid))]];
        v3 = [v3,[TI(seamnode_overlap(:,9)+(j*Ngrid));TI(seamnode_overlap(:,10)+(j*Ngrid))]];
        v2 = [v2,[realnew(tind0+(j*Noverlap));realnew(tind1+(j*Noverlap))]];
    end    
    %% check wether realnew(tind0+(j*Noverlap)) == TI
    
    costAC = CutCost(v1,v2);
    costBC = CutCost(v3,v2);
    costAB = CutCost(v1,v3);
    
%     %% the second method
% 
%     costac = zeros(Nadd,1);
%     costbc = zeros(Nadd,1);
%     costab = zeros(Nadd,1);
%     for i = 0:nbvar-1
%         indexA_0 = seamnode_overlap(:,7) + (i*Ngrid);
%         indexA_1 = seamnode_overlap(:,8) + (i*Ngrid);
%         
%         indexB_0 = seamnode_overlap(:,9) + (i*Ngrid);
%         indexB_1 = seamnode_overlap(:,10) + (i*Ngrid);
%         
%         indexC_0 = tind0 + (i*Noverlap);
%         indexC_1 = tind1 + (i*Noverlap);
%         
%         if variabletype(i+1) == 1
%             tac = abs(TI(indexA_0)-realnew(indexC_0)) + abs(TI(indexA_1)-realnew(indexC_1));
%             tbc = abs(TI(indexB_0)-realnew(indexC_0)) + abs(TI(indexB_1)-realnew(indexC_1));
%             tab = abs(TI(indexA_0)-TI(indexB_0)) + abs(TI(indexA_1)-TI(indexB_1));
%         else
%             im_term1_ac = abs(TI(indexA_0)-realnew(indexC_0)); im_term1_ac(im_term1_ac ~= 0) = 1;
%             im_term2_ac = abs(TI(indexA_1)-realnew(indexC_1)); im_term2_ac(im_term2_ac ~= 0) = 1;
%             
%             im_term1_bc = abs(TI(indexB_0)-realnew(indexC_0)); im_term1_bc(im_term1_bc ~= 0) = 1;
%             im_term2_bc = abs(TI(indexB_1)-realnew(indexC_1)); im_term2_bc(im_term2_bc ~= 0) = 1;
%             
%             im_term1_ab = abs(TI(indexA_0)-TI(indexB_0)); im_term1_ab(im_term1_ab ~= 0) = 1;
%             im_term2_ab = abs(TI(indexA_1)-TI(indexB_1)); im_term2_ab(im_term2_ab ~= 0) = 1;
%             
%             tac = im_term1_ac + im_term2_ac;
%             tbc = im_term1_bc + im_term2_bc;
%             tab = im_term1_ab + im_term2_ab;
%         end
%             costac = costac + (w_v(i+1)*tac); 
%             costbc = costbc + (w_v(i+1)*tbc); 
%             costab = costab + (w_v(i+1)*tab); 
%     end
%     
%     costac = costac + eps;
%     costbc = costbc + eps;
%     costab = costab + eps;
%     
%     %% test 
%     vac = costAC - costac;
%     vbc = costBC - costbc;
%     vab = costAB - costab;
%     
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
    if sum(ACRV == costAB) ~= Nadd
        disp('check cost'); % the difference is 4*eps for each unequal element
    end
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
label = reshape(labels(1:Noverlap),ox,oy);
[AC,Ne] = Addconnected2D(label);

figure(21)
imagesc(label)
axis xy equal tight

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
    end
    
%     for nn = 1:Nadd
%         ch = [num2str(check(nn,1)),num2str(check(nn,2)),num2str(check(nn,3))];            
%         switch ch
%             case '000'
%                 % keep the old cut
%                 keep_old = [keep_old;seamnode_overlap(nn,:)];
%             case '110'
%                 % make a new cut at the same place
%                 % remove the old one and keep the new one
%  %               new_i = seamnode_overlap(nn,:);
%                 changeline1 = find(seam_new(:,5)==seamnode_overlap(nn,6) & seam_new(:,6)==seamnode_overlap(nn,5));
%                 changeline = find(new_ind0(:)==tind1(nn) & new_ind1(:)==tind0(nn));
%                 if changeline1 ~= changeline
%                     disp('unequal line 147');
%                 end
%                 seam_new(changeline,4) = costBC(nn);
%                 seam_new(changeline,8) = seamnode_overlap(nn,9);
%             case '101'
%                 % make a new cut at the same place
%                 % remove the old one and keep the new one
%                 changeline1 = find(seam_new(:,5)==seamnode_overlap(nn,5) & seam_new(:,6)==seamnode_overlap(nn,6));
%                 changeline = find(new_ind0(:)==tind0(nn) & new_ind1(:)==tind1(nn));
%                 if changeline1 ~= changeline
%                     disp('unequal line 156');
%                 end
%                 seam_new(changeline,4) = costAC(nn);
%                 seam_new(changeline,8) = seamnode_overlap(nn,8);
%             case '111'
%                 % remove the old cut
%             case '001' % AC
% %                 troubleline1 = find(seam_new(:,5)==seamnode_overlap(nn,5) & seam_new(:,6)==seamnode_overlap(nn,6));
% %                 troubleline = find(new_ind0(:)==tind0(nn) & new_ind1(:)==tind1(nn));
% %                 c = [costAB(nn),costAC(nn),costBC(nn)];
% %                 seamnew_trouble = seam_new(troubleline,:);
% %                 seamold = seamnode_overlap(nn,:);
%                 % AB + BC < AC
%                 % make a new cut at the same place
%                 % remove the old one and keep the new one
%                 changeline1 = find(seam_new(:,5)==seamnode_overlap(nn,5) & seam_new(:,6)==seamnode_overlap(nn,6));
%                 changeline = find(new_ind0(:)==tind0(nn) & new_ind1(:)==tind1(nn));
%                 if changeline1 ~= changeline
%                     disp('unequal line 156');
%                 end
%                 seam_new(changeline,4) = costAC(nn);
%                 seam_new(changeline,8) = seamnode_overlap(nn,8);
%                 disp(ch);
%             case '010' %BC
% %                 troubleline = find(new_ind0(:)==tind1(nn) & new_ind1(:)==tind0(nn));
% %                 troubleline1 = find(seam_new(:,5)==seamnode_overlap(nn,6) & seam_new(:,6)==seamnode_overlap(nn,5));
% %                 c = [costAB(nn),costAC(nn),costBC(nn)];
% %                 seamnew_trouble = seam_new(troubleline,:);
% %                 seamold = seamnode_overlap(nn,:);
%                 % AB + AC < BC
%                 % make a new cut at the same place
%                 % remove the old one and keep the new one
%  %               new_i = seamnode_overlap(nn,:);
%                 changeline1 = find(seam_new(:,5)==seamnode_overlap(nn,6) & seam_new(:,6)==seamnode_overlap(nn,5));
%                 changeline = find(new_ind0(:)==tind1(nn) & new_ind1(:)==tind0(nn));
%                 if changeline1 ~= changeline
%                     disp('unequal line 147');
%                 end
%                 seam_new(changeline,4) = costBC(nn);
%                 seam_new(changeline,8) = seamnode_overlap(nn,9);
%                 disp(ch);
%             case '100' 
% %                 troubleline1 = find(new_ind0(:)==tind1(nn) & new_ind1(:)==tind0(nn));
% %                 troubleline2 = find(new_ind0(:)==tind0(nn) & new_ind1(:)==tind1(nn));
% %                 troubleline = [troubleline1,troubleline2];
% %                 c = [costAB(nn),costAC(nn),costBC(nn)];
% %                 seamnew_trouble = seam_new(troubleline,:);
% %                 seamold = seamnode_overlap(nn,:);
%                 keep_old = [keep_old;seamnode_overlap(nn,:)];
%                 disp(ch);
%                 % AC + BC < AB
%             otherwise
%                 disp('unexpected cut!'); 
%                 disp(ch);
%                 troublelineAC = find(new_ind0(:)==tind0(nn) & new_ind1(:)==tind1(nn));
%                 troublelineBC = find(new_ind0(:)==tind1(nn) & new_ind1(:)==tind0(nn));
%                 troubleline = [troublelineAC;troublelineBC];
%                 c = [costAB(nn),costAC(nn),costBC(nn)];
%                 if size(troubleline,1) > 0
%                     seamnew_trouble = seam_new(troubleline,:);
%                 end
%                 seamold = seamnode_overlap(nn,:);
%                 subtract = costAB(nn)- seamnode_overlap(nn,4);
%         end
%     end

end

%newseamnode = seamnode;
newseamnode = [seamnode;seam_new;keep_old1];

costmap = zeros(SGdim(1:2));
costmap(newseamnode(:,5)) = 0.5 * newseamnode(:,4);
costmap(newseamnode(:,6)) = 0.5 * newseamnode(:,4);

figure(22); clf;
imagesc(costmap); axis xy equal tight;

end