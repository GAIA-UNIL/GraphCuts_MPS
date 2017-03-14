function [modiSG,labelmap,strucmap,ch] = modi_con_3D(rr,TI,condSG,seamnode,costmap,labelmap,strucmap,filename)

global nbreplicates stop maxloop
global indexmap SGdim halfpr

meanerr = sum(seamnode(:,4))/size(seamnode,1);
ch = meanerr;
l = 0; r = 1.3;
sigb = ones(SGdim); sigb(2:SGdim(1)-1,2:SGdim(2)-1,2:SGdim(3)-1)=0;
start = max(strucmap(:));

while (l < maxloop && meanerr>stop)
    %% sort the costmap
    costv = costmap(costmap>eps);
    costv = sort(costv(:),'descend');
    threadhold = quantile(costv(:), 0.8);
    mod = zeros(size(costmap));
    mod(costmap(:) >= threadhold) = 1;
    
    [pat, nri] = bwlabeln(mod,6);
    if nri >1
        nbnr = histc(pat(:),[1:nri]);
        [nbnri,i] = sort(nbnr(:),'descend');
        
        if nri >= nbreplicates
            selec = i(randi(nbreplicates));
        else
            selec = i(randi(nri));
        end
    else
        selec = 1;
    end
    
    pati = zeros(SGdim);
    pati(pat == selec) = 1;
    patI = find(pati == 1);
    [patdimx,patdimy,patdimz] = ind2sub(SGdim,patI);
    
    mindimx = min(patdimx);    maxdimx = max(patdimx);
    mindimy = min(patdimy);    maxdimy = max(patdimy);
    mindimz = min(patdimz);    maxdimz = max(patdimz);
    
    centerx = floor((mindimx+maxdimx)/2);
    centery = floor((mindimy+maxdimy)/2);
    centerz = floor((mindimz+maxdimz)/2);
    
    presx = maxdimx - mindimx;
    presy = maxdimy - mindimy;
    presz = maxdimz - mindimz;
    
    halfx = min(ceil(presx * r /2),floor(SGdim(1)*0.5)-1);
    halfx = max(halfx,halfpr(1));
    halfy = min(ceil(presy * r /2),floor(SGdim(2)*0.5)-1);
    halfy = max(halfy,halfpr(2));
    halfz = min(ceil(presz * r /2),floor(SGdim(3)*0.5)-1);
    halfz = max(halfz,halfpr(3));
    
    SGpickx = max(centerx - halfx, 1);
    SGpicky = max(centery - halfy, 1);
    SGpickz = max(centerz - halfz, 1);
    
    patchmaxx = min(centerx + halfx,SGdim(1));
    patchmaxy = min(centery + halfy,SGdim(2));
    patchmaxz = min(centerz + halfz,SGdim(3));
    
    patchsizex = patchmaxx - SGpickx + 1;
    patchsizey = patchmaxy - SGpicky + 1;
    patchsizez = patchmaxz - SGpickz + 1;
    
    realold = condSG(SGpickx:patchmaxx,SGpicky:patchmaxy,SGpickz:patchmaxz);
    labelold = labelmap(SGpickx:patchmaxx,SGpicky:patchmaxy,SGpickz:patchmaxz);
    sigbi = sigb(SGpickx:patchmaxx,SGpicky:patchmaxy,SGpickz:patchmaxz);
    strucmapi = strucmap(SGpickx:patchmaxx,SGpicky:patchmaxy,SGpickz:patchmaxz);
    
    %% select a new patch
%    TIs = TI(1:TIdim(1)-patchsizex+1,1:TIdim(2)-patchsizey+1,1:TIdim(3)-patchsizez+1);
    dis = dis3D(TI,realold);
    [v,i]=sort(dis(:));
    if size(i,1) >= nbreplicates
       sel = randi(nbreplicates);
    else
        sel = randi(size(i,1));
    end
    ind = i(sel);
    [TIpickx,TIpicky,TIpickz]=ind2sub(size(dis),ind);
    realnew = TI(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1,TIpickz:TIpickz+patchsizez-1);
    labelnew = indexmap(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1,TIpickz:TIpickz+patchsizez-1);

    errori = abs(realold - realnew);
    
    %% define terminal
    Tb = ones(patchsizex,patchsizey,patchsizez);
    Tb(2:patchsizex-1,2:patchsizey-1,2:patchsizez-1)=0;
    Td = find(sigbi == 1); 
    T1 = find(Tb == 1);
    T1 = setdiff(T1,Td);
    
    pro = 0.8;
    threadpro = quantile(errori(:),pro);
    c = errori(errori ~= 0);
    if threadpro < mean(c(:))
        threadpro = mean(c(:));
    end
    paterr = (errori>=threadpro);
    paterr = double(paterr);
    paterr(T1) = 0;
    [patl,patn] = bwlabeln(paterr,6);
    nbpatn = histc(patl(:),[1:patn]);
    [nbpatni,i] = sort(nbpatn(:), 'descend');
    if patn >= nbreplicates
        selt = i(randi(nbreplicates));
    else if patn<nbreplicates && patn>=1
            selt = i(randi(patn));
        else
            l = l + 1;
            continue;
        end
    end
    T2 = find(patl == selt);
    T2 = setdiff(T2,T1);
    
    %% make a new cut
    [oldseamnode,newseamnode,label,preseam] = AddseamCut3D(errori,seamnode,TI, realnew,labelold, labelnew,SGpickx,SGpicky,SGpickz,T1,T2);
    if sum(newseamnode(:,4)) <= sum(oldseamnode(:,4))
        % make a new cut
        seamnode = newseamnode;
        realprop = realnew;
        realprop(label == 0) = realold(label==0);
        condSG(SGpickx:patchmaxx,SGpicky:patchmaxy,SGpickz:patchmaxz) = realprop;
        
        labelprop = labelnew;
        labelprop(label == 0) = labelold(label == 0);
        labelmap(SGpickx:patchmaxx,SGpicky:patchmaxy,SGpickz:patchmaxz)=labelprop;
        
        start = start + 1;
        strucmapi(label==1) = start;
        strucmap(SGpickx:patchmaxx,SGpicky:patchmaxy,SGpickz:patchmaxz) = strucmapi;
        
        costmap = zeros(SGdim(1:3));
        costmap(seamnode(:,5)) = 0.5 * seamnode(:,4);
        costmap(seamnode(:,6)) = 0.5 * seamnode(:,4);
        
        meanerr = sum(seamnode(:,4))/size(seamnode,1);
        ch = [ch;meanerr];
    else
        seamnode = oldseamnode;
    end
    l = l+1;
end

if size(ch,1) > 1
    modiSG = condSG;
    condnamei = ['modi',filename,num2str(rr),'.SGEMS'];
    WriteGrid(modiSG,condnamei,'modi');
    WriteGrid(strucmap,['modi-struct-',num2str(rr),'.SGEMS'],'modi');
else
    modiSG = nan;
end
figure(1); clf;
plot([1:size(ch,1)],ch);
% saveas(gcf,'err-change.fig');

end
