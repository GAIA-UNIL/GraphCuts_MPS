function [modiSG,meancost,seamnode] = modi_2D_con(rr,condSG,TI,SGname,seamnode,Cond,costmap,errorco,labelmap)
% the function to modify the initial realization according to the costmap
global nbreplicates
global showfig
global indexmap
global SGdim
global w
global nbvar
global variabletype
global maxloop
global stop
global patchrange
global r p th

stopping = stop * max(costmap(:));
meanseam = mean(seamnode(:,4));
sigb = ones(SGdim(1),SGdim(2)); sigb(2:SGdim(1)-1,2:SGdim(2)-1) = 0;
meancost = [];
modiSG = condSG;

ll = 0; step = 1;
while ll<maxloop && meanseam>stopping
    %% sort the costmap
    costv = costmap(costmap>eps);
    costv = sort(costv(:),'descend');
    threadhold = quantile(costv(:), 0.8);
    
    mod = zeros(size(costmap));
    mod(costmap(:) >= threadhold) = 1;
    
    [pat, nri] = bwlabel(mod,4);
    if nri > 1
        nbnr = histc(pat(:),[1:nri]);
        [nbnri,i] = sort(nbnr(:), 'descend');
        
        if nbreplicates>1 && nri>=nbreplicates
            selec = i(randi(nbreplicates));
        else
            selec = 1;
        end
    else
        selec = 1;
    end
    
    pati = zeros(size(costmap));
    pati(pat == selec) = 1;
    [patdimx, patdimy] = find(pati == 1);
    patin = sub2ind(size(pati), patdimx,patdimy);
    mindimx = min(patdimx);
    maxdimx = max(patdimx);
    mindimy = min(patdimy);
    maxdimy = max(patdimy);
    
    centerx = floor((mindimx+maxdimx)/2);
    centery = floor((mindimy+maxdimy)/2);
    
    patchsizex = maxdimx - mindimx;
    patchsizey = maxdimy - mindimy;
%     halfdimx = min(ceil(patchsizex * r /2),patchrange(1,2));
%     halfdimy = min(ceil(patchsizey * r /2),patchrange(2,2));
    halfdimx = min(ceil(patchsizex * r /2),floor(SGdim(1)*0.5)-1);
    halfdimy = min(ceil(patchsizey * r /2),floor(SGdim(2)*0.5)-1);
    halfdimx = max(halfdimx,patchrange(1,1));
    halfdimy = max(halfdimy,patchrange(2,1));
    
    SGpickx = max(centerx - halfdimx, 1);
    SGpicky = max(centery - halfdimy, 1);
    patchmaxx = min(centerx + halfdimx,SGdim(1));
    patchmaxy = min(centery + halfdimy,SGdim(2));
    patchsizex = patchmaxx - SGpickx + 1;
    patchsizey = patchmaxy - SGpicky + 1;
    
    realold = modiSG(SGpickx:patchmaxx,SGpicky:patchmaxy,:);
    labelold = labelmap(SGpickx:patchmaxx,SGpicky:patchmaxy);
    sigforcondi = Cond(SGpickx:patchmaxx,SGpicky:patchmaxy,:);
    
    sigbi = sigb(SGpickx:patchmaxx,SGpicky:patchmaxy);
    delb = find(sigbi == 1);
    % the boundary of SG should be deleted from the terminal
    
    %% select a new patch
%     uncon_dis = distance2D(TI,realold);
%     con_dis = distance2D(TI,sigforcondi);
%     filtermap = ((1-w)* uncon_dis) + (w*con_dis);
    
    filtermap = dis_con2D_v2(TI,realold,sigforcondi,w);
    [v,i] = sort(filtermap(:));
    if nbreplicates > 1
        sel = i(randi(nbreplicates));
    else
        sel = i(1);
    end
    [fx,fy] = size(filtermap);
    [TIpickx,TIpicky] = ind2sub([fx,fy],sel);
    realnew = TI(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1,:);
    labelnew = indexmap(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1);
    
    errori = multidis(realnew,realold);
    
    %% Define terminal
    Tbd = ones(patchsizex,patchsizey);
    Tbd(2:patchsizex-1,2:patchsizey-1) = 0;
    Tb = find(Tbd == 1);
    
    pro = 0.8;
    threadpro = quantile(errori(:),pro);
    if threadpro < mean(errori(:))
        threadpro = mean(errori(:));
    end
    paterr = (errori>=threadpro);
    paterr = double(paterr);
    %paterr(T1) = 0;
    [patl,patn] = bwlabel(paterr,4);
    nbpatn = histc(patl(:),[1:patn]);
    [nbpatni,i] = sort(nbpatn(:), 'descend');
    if patn >= nbreplicates
        selt = randi(nbreplicates);
    else if patn<nbreplicates && patn>1
            selt = randi(patn);
        else if patn == 1
                selt = 1;
            else
                ll = ll + 1;
                continue;
            end
        end
    end
    patv = i(selt);
    T2 = find(patl == patv);
    
    %% Keep the conditioning data
    c_It = find(isfinite(sigforcondi));
    %% detect conditioning
    checkIc = abs(realold(c_It) - realnew(c_It));
    delIc = find(checkIc <= th);
    c_It(delIc) = [];
    [Icx_i,Icy_i,Ind] = ind2sub(size(realnew),c_It);
    Ic_i = sub2ind([patchsizex,patchsizey],Icx_i,Icy_i);
    Ic_i = unique(Ic_i);
    nbTc = size(Ic_i,1);    
    T1_c = [];
    
    if nbTc > 0
        for uu = 1:nbTc
            rangmin = errori(Ic_i(uu))*(1-p);
            rangmax = errori(Ic_i(uu))*(1+p);
            sigerr = (errori>=rangmin & errori<=rangmax);
            fsig = bwlabel(sigerr,4);
            T1_enc = find(fsig == fsig(Ic_i(uu)));
            T1_c = unique([T1_c; T1_enc]);
        end
    end
    T1 = union(Tb,T1_c);
    delb = setdiff(delb,Ic_i);
    T1 = setdiff(T1,delb);
    
    % sink is defined according to costmap
    pro = 0.8;
%     non_e = errori;
%     non_e(non_e == 0) = [];
%     threadpro1 = quantile((non_e(:)),pro);
    
    threadpro = quantile(errori(:),pro);
    if threadpro < mean(errori(:))
        threadpro = mean(errori(:));
    end
    paterr = (errori>=threadpro);
    paterr = double(paterr);
    paterr(T1) = 0;
    [patl,patn] = bwlabel(paterr,4);
    nbpatn = histc(patl(:),[1:patn]);
    [nbpatni,i] = sort(nbpatn(:), 'descend');
    if patn >= nbreplicates
        selt = randi(nbreplicates);
    else if patn<nbreplicates && patn>1
            selt = randi(patn);
        else if patn == 1
                selt = 1;
            else
                ll = ll + 1;
                continue;
            end
        end
    end
    patv = i(selt);
    T2 = find(patl == patv);
    T2 = setdiff(T2,T1);
    Np = patchsizex * patchsizey;
    
    %% Scan the seam nodes and make a new cut
    [oldseamnode,newseamnode,label,preseam] = AddseamCut2D(errori,seamnode,TI, realnew,labelold, labelnew,SGpickx,SGpicky,T1,T2);
    nb_oldseam = size(oldseamnode,1);
    nb_newseam = size(newseamnode,1);
    
    eva_old_cost = sum(oldseamnode(:,4)) / nb_oldseam;
    eva_new_cost = sum(newseamnode(:,4)) / nb_newseam;
    
    %% Update
    if eva_new_cost <= eva_old_cost
        Np = patchsizex * patchsizey;
        seamnode = newseamnode;
        realprop = realnew;
        changeind = find(label==0);
        for k = 0:nbvar-1
            change = changeind + k*Np;
            realprop(change) = realold(change);
        end         
        modiSG(SGpickx:patchmaxx,SGpicky:patchmaxy,:) = realprop;
        
        labelprop = labelnew;
        labelprop(changeind) = labelold(changeind);
        labelmap(SGpickx:patchmaxx,SGpicky:patchmaxy) = labelprop;
        
        costmap = zeros(SGdim(1:2));
        costmap(seamnode(:,5)) = 0.5 * seamnode(:,4);
        costmap(seamnode(:,6)) = 0.5 * seamnode(:,4);
        
        step = step + 1;
        if showfig == 1
            Td = zeros(patchsizex,patchsizey);
            Td(T1) = 1; Td(T2) = 2;
            showfigures(label,Td,preseam,costmap,modiSG,rr,step,'_modi');
        end
        
        meanseam = eva_new_cost;
        meancost = [meancost;eva_new_cost];
    else
        seamnode = oldseamnode;
    end
    
    ll = ll + 1;
 end
end