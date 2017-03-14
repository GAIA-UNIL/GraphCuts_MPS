function [modiSG,meancost,seamnode] = modi_2D(rr,initSG,TI,outputname,seamnode,costmap,labelmap)
% the function to modify the initial realization according to the costmap
global nbreplicates
global showfig
global indexmap
global SGdim
global w_v
global nbvar
global variabletype
global maxloop
global stop
global patchrange
global r

stopping = stop * max(costmap(:));
meanseam = mean(seamnode(:,4));
sigb = ones(SGdim(1),SGdim(2)); sigb(2:SGdim(1)-1,2:SGdim(2)-1) = 0;
meancost = [];

ll = 0; step = 1;
while ll<maxloop && meanseam>stopping
%     patchsizex = patchrange(1,1) + randi(patchrange(1,2)-patchrange(1,1));
%     patchsizey = patchrange(2,1) + randi(patchrange(2,2)-patchrange(2,1));
%     halfx = round(patchsizex);
%     halfy = round(patchsizey);
%     
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
    
    realold = initSG(SGpickx:patchmaxx,SGpicky:patchmaxy,:);
    labelold = labelmap(SGpickx:patchmaxx,SGpicky:patchmaxy);
    
    sigbi = sigb(SGpickx:patchmaxx,SGpicky:patchmaxy);
    delb = find(sigbi == 1);
    % the boundary of SG should be deleted from the terminal
    
    %% select a new patch
     filtermap = distance2D(TI,realold);
%    [filtermap,No] = distance2D_v1(TI,realold);
    [v,i] = sort(filtermap(:));
    if nbreplicates > size(i,1)
        sel = i(randi(nbreplicates));
    else
        sel = i(randi(size(i,1)));
    end
    [fx,fy] = size(filtermap);
    [TIpickx,TIpicky] = ind2sub([fx,fy],sel);
    realnew = TI(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1,:);
    labelnew = indexmap(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1);
    
    errori = zeros(patchsizex,patchsizey);
    for k = 1:nbvar
        if variabletype(k) == 1
            errori = errori + w_v(k)*abs(realnew(:,:,k)-realold(:,:,k));
        else
            t_ei = abs(realnew(:,:,k)-realold(:,:,k));
            t_ei(t_ei ~= 0) = 1;
            errori = errori + w_v(k)*t_ei;
        end
    end
    
    %% Define terminal
    Tb = ones(patchsizex,patchsizey);
    Tb(2:patchsizex-1,2:patchsizey-1) = 0;
    T1 = find(Tb == 1);
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
%     [lineno,flowold,flow,labels,Nadd,check,preseam] = Addseam2D(errori,seamlist,TI,realnew,SGpickx,SGpicky,T1,T2);
%     label = reshape(labels(1:Np),patchsizex,patchsizey);
%     [AC,Ne] = Addconnected2D(label);
%     newseamnode = newseam(AC,ii,jj,1,labelrealnew,labelrealold, errori,TI);
%     [seamlist,costmap] = updateseam2D(seamlist,newseamnode,costmap,Nadd,lineno,check);

%    [oldseamnode,newseamnode,label,preseam] = AddseamCut2D_v2(errori,seamnode,TI, realnew,labelold, labelnew,SGpickx,SGpicky,T1,T2);
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
        initSG(SGpickx:patchmaxx,SGpicky:patchmaxy,:) = realprop;
        
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
            showfigures(label,Td,preseam,costmap,initSG,rr,step,'_modi');
        end
        
        meanseam = eva_new_cost;
        meancost = [meancost;eva_new_cost];
    else
        seamnode = oldseamnode;
    end
    
    ll = ll + 1;

end

modiSG = initSG;
modinamei = ['modi',outputname,num2str(rr),'.SGEMS'];
WriteGrid(modiSG,modinamei,'modi');