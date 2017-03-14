function [lineno, flowold, flownew, meanolderr, labels,Nadd,Ng, check, inpatchlist0, inpatchlist1] = Addseamnodes_3D(seamlist,adjacent,startx,starty,startz, simoldstay,errori,T1,T2,SGdim)

[simsizex,simsizey,simsizez] = size(simoldstay);
Ng = simsizex*simsizey*simsizez;
lineno1 = find(seamlist(:,1)>startx & seamlist(:,1)<startx+simsizex-1);
lineno2 = find(seamlist(:,2)>starty & seamlist(:,2)<starty+simsizey-1);
lineno3 = find(seamlist(:,3)>startz & seamlist(:,3)<startz+simsizez-1);
linenoxy = intersect(lineno1,lineno2);
lineno = intersect(linenoxy,lineno3);

if size(lineno,1)*size(lineno,2) == 0
    N = Ng;
    E = edges6connected3D(simoldstay);
    V = errori(E(:,1)) + errori(E(:,2)) + eps;
    Nadd = 0;
    flowold = 0;
    meanolderr = 0;
    inpatchlist0 = [];
    inpatchlist1 = [];
else
    Nadd = size(lineno,1);
    N = Ng + Nadd;
    addno = [Ng+1:N]';
    size2 = size(T2,1);
    if size2 > 1
        Ttemp = T2(randi(size2));
    else
        Ttemp = T2(1);
    end
    
    Edif = nan(Nadd*6, 2);
    Vdif = nan(Nadd*6, 1);
    
    Esame = edges6connected3D(simoldstay);
    
    %% case 1: edges between seam node and sink
    Edif(1:Nadd,1) = addno;
    Edif(1:Nadd,2) = Ttemp;
    Edif(Nadd+1:2*Nadd, 1) = Ttemp;
    Edif(Nadd+1:2*Nadd, 2) = addno;
    
    %% cost of edge case 1
    Vdif(1:Nadd) = seamlist(lineno,4);
    Vdif(Nadd+1:2*Nadd) = seamlist(lineno,4);
    
    %% case 2: edges between seam node and label 0 adjacent node
    list0 = adjacent(lineno,1);
    [sx,sy,sz] = ind2sub([SGdim(1),SGdim(2),SGdim(3)],list0);
    inpatchlist0 = sub2ind([simsizex,simsizey,simsizez],sx-startx+1,sy-starty+1,sz-startz+1);
    
    Edif(2*Nadd+1:3*Nadd,1) = inpatchlist0;
    Edif(2*Nadd+1:3*Nadd,2) = addno;
    Edif(3*Nadd+1:4*Nadd,1) = addno;
    Edif(3*Nadd+1:4*Nadd,2) = inpatchlist0;
    
    
    %% case 3: edges between seam node and label 1 adjacent node
    list1 = adjacent(lineno,2);
    [sx,sy,sz] = ind2sub([SGdim(1),SGdim(2),SGdim(3)],list1);
    inpatchlist1 = sub2ind([simsizex,simsizey,simsizez],sx-startx+1,sy-starty+1,sz-startz+1);
    
    Edif(4*Nadd+1:5*Nadd,1) = addno;
    Edif(4*Nadd+1:5*Nadd,2) = inpatchlist1;
    Edif(5*Nadd+1:6*Nadd,1) = inpatchlist1;
    Edif(5*Nadd+1:6*Nadd,2) = addno;
    
    %% cost of edge case 2
    err1 = abs(adjacent(lineno,3) - simoldstay(inpatchlist1));
    err2 = abs(simoldstay(inpatchlist0) - adjacent(lineno,4));
    Vdif(2*Nadd+1:3*Nadd) = err1 + err2 + eps;
    Vdif(3*Nadd+1:4*Nadd) = Vdif(2*Nadd+1:3*Nadd);
    
    %% cost of edge case 3
    err1_2 = abs(simoldstay(inpatchlist0) - adjacent(lineno,6));
    err2_2 = abs(adjacent(lineno,5) - simoldstay(inpatchlist1));
    Vdif(4*Nadd+1:5*Nadd) = err1_2 + err2_2 + eps;
    Vdif(5*Nadd+1:6*Nadd) = err1_2 + err2_2 + eps;
    
    %% remove edges in Esame
    remove = nan(2*Nadd,2);
    remove(:,1) = [inpatchlist0;inpatchlist1];
    remove(:,2) = [inpatchlist1;inpatchlist0];
    for uu = 1:2*Nadd
        removeline = find(Esame(:,1)==remove(uu,1) & Esame(:,2)==remove(uu,2));
        Esame(removeline,:) = [];
    end
    Vsame = errori(Esame(:,1)) + errori(Esame(:,2)) + eps;
    
    E = [Esame;Edif];
    V = [Vsame;Vdif];
    
    flowoldlist = seamlist(lineno,4);
    flowold = sum(flowoldlist(:));
    meanolderr = flowold/Nadd;
    
end

size1 = size(T1,1);size2=size(T2,1);
A = sparse(E(:,1),E(:,2),V,N,N);
T = sparse([T1;T2],[ones(size1,1);ones(size2,1)*2],ones(size1+size2,1)*9e9, N ,2);
[flow, labels] = maxflow(A,T);

flownew = flow;
%% check the value of the added seam nodes
check = nan(Nadd,3); 
if Nadd > 0
    check(:,1) = labels(Ng+1:N);
    check(:,2) = labels(inpatchlist0);
    check(:,3) = labels(inpatchlist1);    
end


end