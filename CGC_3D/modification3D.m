function [ SGi,seamlist, adjacent, nbvalid, meanseam ] = modification3D(SG,TI,sr,halfrange,nbreplicate,seamlist,adjacent,costmap,nbloops, maxloop, stop, ep, showfig)
%MODIFICATION3D Summary of this function goes here
%   Detailed explanation goes here
[SGdimx,SGdimy,SGdimz] = size(SG);
NSG = SGdimx*SGdimy*SGdimz;
SGi = nan(NSG,nbloops);
meanseam = nan(nbloops,1);
maxerr = max(seamlist(:,4));
nbvalid = 0;
nbcount = 1;
sigb = ones(SGdimx,SGdimy,SGdimz); sigb(2:SGdimx-1,2:SGdimy-1,2:SGdimz-1)=0;
while maxerr>stop || £¨nbvalid<nbloops && nbcount<=maxloop£©
    
    %% sort the costmap
    costv = costmap(costmap>0);
    costv = sort(costv(:),'descend');
    threadhold = quantile(costv(:),0.7);
    
    mod = zeros(SGdimx,SGdimy,SGdimz);
    mod(costmap(:)>=threadhold) = 1;
    figure(11);
    ViewGrid(mod); view(30,30);
    [pat,nri] = bwlabel(mod,4);
    
    if nri>1
        nbnr = histc(pat(:),[1:nri]);
        [nbnri,i] = sort(nbnr(:),'descend');
        
        if nbreplicate>1 && nri>nbreplicate
            selec = i(randi(nbreplicate));
        else
            selec = 1;
        end
    else
        selec = 1;
    end
    patind = find(pat == selec);
    [patx,paty,patz]=sub2ind([SGdimx,SGdimy,SGdimz], patind);
    minx = min(patx);
    miny = min(paty);
    minz = min(patz);
    maxx = max(patx);
    maxy = max(paty);
    maxz = max(patz);
    
    centerx = floor((minx+maxx)/2);
    centery = floor((miny+maxy)/2);
    centerz = floor((minz+maxz)/2);
    
    selhalfx = ceil((maxx - minx)*sr/2);
    selhalfy = ceil((maxy - miny)*sr/2);
    selhalfz = ceil((maxz - minz)*sr/2);
    
    selhalfx = max(halfrange(1,1),selhalfx);
    selhalfx = min(selhalfx,halfrange(1,2));
    selhalfy = max(selhalfy,halfrange(2,1));
    selhalfy = min(selhalfy,halfrange(2,2));
    selhalfz = max(selhalfz,halfrange(3,1));
    selhalfz = min(selhalfz,halfrange(3,2));
    
    SGpickx = max(centerx-selhalfx,1);
    SGpickmaxx = min(centerx+selhalfx,SGdimx);
    SGpicky = max(centery-selhalfy,1);
    SGpickmaxy = min(centery+selhalfy,SGdimy);
    SGpickz = max(centerz-selhalfz,1);
    SGpickmaxz = min(centerz+selhalfz,SGdimz);
    
    patchsizex = SGpatchmaxx-SGpickx+1;
    patchsizey = SGpatchmaxy-SGpicky+1;
    patchsizez = SGpatchmaxz-SGpickz+1;
    
    Noverlap = patchsizex*patchsizey*patchsizez;
    %% old patch
    realold = SG(SGpickx:SGpickmaxx,SGpicky:SGpickmaxy,SGpickz:SGpickmaxz);
    sigi = sigb(SGpickx:SGpickmaxx,SGpicky:SGpickmaxy,SGpickz:SGpickmaxz);
    
    %% search training image to find a new patch
    filtermap = distance3D(realold,TI);
    selt = randi(nbreplicate);
    [v,i] = sort(filtermap(:));
    seltind = i(selt);
    [TIpickx,TIpicky,TIpickz] = ind2sub(size(filtermap),seltind);
    realnew = TI(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1,TIpickz:TIpickz+patchsizez-1);
    
    %% define terminal
    errori = abs(realold - realnew);
    threadp = quantitle(errori(:),ep);
    paterr = (errori>threadp);   paterr = double(paterr);
    [paterrv,paterrn] = bwlabel(paterr,4);
    if paterrn >= nbreplicate
        sele = randi(nbreplicate);
    else if paterrn>1 && paterrn<nbreplicate
            sele = randi(nbreplicate);
        else if paterrn == 1
                sele = 1;
            else
                nbcount = nbcount + 1;
                continue;
            end
        end
    end
        
    nbpatn = histc(paterrv(:),[1:paterrn]);
    [nbpatni,paterri] = sort(nbpatn(:),'descend');
    selpaterrv = paterri(sele);
    T2e = find(paterrv == selpaterrv);    
    
    delb = find(sigbi == 1);
    sigp = ones(patchsizex,patchsizey,patchsizez); sigp(2:patchsizex-1,2:patchsizey-1,2:patchsizez-1)=0;
    T1_b = find(sigp == 1);
    
    T1 = setdiff(TI_b,delb);
    T2 = setdiff(T2e,T1_b);
    
    %% search seamlist to find old cuts
    [lineno, flowold, flownew, meanolderr, labels,Nadd,Ng, check, inpatchlist0, inpatchlist1] = Addseamnodes_3D(seamlist,adjacent,startx,starty,startz, simoldstay,errori,T1,T2,SGdim);
    label = reshape(labels(1:Noverlap));
    if Nadd >0
        listold = seamlist(:,4);
        evaold = sum(listold(:))/Nadd;
        maxoldcost = max(listold(:));
    else
        evaold = 0;
        maxoldcost = 0;
    end
    %% update seamlist
	[seamlistnew,adjacentnew] = updateseam3D(lineno,seamlist, adjacent, inpatchlist0, inpatchlist1, check, label,ii,jj,kk,realold,realnew, Nadd,SGdim);
    Naddnew = size(seamlistnew,1);
    if Naddnew > 0
        seamvalue = seamlistnew(:,4);
        evanew = sum(seamvalue(:))/Naddnew;
        maxnewcost = max(seamvalue(:));
    else
        evanew = 0;
        maxnewcost = 0;
    end
    
    %% rejection strategy
    if evanew <= evaold %accept
        nbvalid = nbvalid + 1;
        nbcount = nbcount + 1;
        %% update seam node
        costmap(adjacent(lineno,1:2)) = 0;
        costmap(adjacentnew(:,1)) = abs(adjacentnew(:,3)-adjacentnew(:,5));
        costmap(adjacentnew(:,2)) = abs(adjacentnew(:,4)-adjacentnew(:,6));
        
        seamlist(lineno,:) = [];
        adjacent(lineno,:) = [];
        seamlist = [seamlistnew;seamlist];
        adjacent = [adjacentnew;adjacent];
        maxerr = max(seamlist(:,4));
        meanseam(nbvalid) = sum(seamlist(:,4))/size(seamlist,1);
        
        %% update simulation grid
        realprop = realnew;
        realprop(label==0) = realold(label == 0);
        SG(SGpickx:SGpickmaxx,SGpicky:SGpickmaxy,SGpickz:SGpickmaxz) = realprop;
        SGprop = reshape(SG,NSG,1);
        SGi(:,nbvalid) = SGprop;
        
        %% show figures
        if showfig == 1
            figure(2);
            ViewGrid(SG);   axis equal tight;   view(30,30);
               
            figure(3)
            ViewPts(seamlist,1,20,'filled');    axis equal tight;
            view(30,30);    title('costmap');
        
            figure(4);  clf
            subplot(2,2,1);
            ViewGrid(realold);
            axis equal tight;   view(30,30);    title('old patch')
            subplot(2,2,2);
            ViewGrid(realnew);
            axis equal tight;   view(30,30);    title('new patch')
            subplot(2,2,3);
            ViewGrid(label);    
            axis equal tight;   view(30,30);    title('cut')
            subplot(2,2,4);
            ViewGrid(realprop);
            axis equal tight;   view(30,30);    title('stitched patch')
        end
    else %reject
        nbcount = nbcount + 1;
        continue;
    end

end
end

