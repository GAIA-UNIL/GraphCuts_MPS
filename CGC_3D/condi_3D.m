function [condSG,seamnode,costmap,labelmap,strucmap,meancoerr] = condi_3D(rr,TI,initSG,conditioning, condierr, filename,seamnode,labelmap,costmap,strucmap)

global indexmap SGdim w_c halfpr nbreplicates TIdim

p=1.3;
indcond = sub2ind(SGdim,conditioning(:,1),conditioning(:,2),conditioning(:,3));
cond = nan(SGdim); cond(indcond)=conditioning(:,4);
sigb = zeros(SGdim); 
sigb(SGdim(1),:,:) = 1; sigb(:,SGdim(2),:) = 1;
a = find(condierr ~= 0);
nbch = size(a,1);
ns = max(strucmap(:));
meancoerr = [];

while nbch > 0
    newm = sum(condierr(:));
    meancoerr = [meancoerr;newm];
    %% choose a un-match conditioning data to modify
    if nbch > 1
        pe = randi(nbch);
        centerI = a(pe);
    else
        centerI = a;
    end
    [centerx,centery,centerz] = ind2sub(SGdim,centerI);
    minx = max(1,centerx-halfpr(1)); maxx = min(centerx+halfpr(1),SGdim(1));
    miny = max(1,centery-halfpr(2)); maxy = min(centery+halfpr(2),SGdim(2));
    minz = max(1,centerz-halfpr(3)); maxz = min(centerz+halfpr(3),SGdim(3));
    
    patchsizex = maxx - minx + 1;
    patchsizey = maxy - miny + 1;
    patchsizez = maxz - minz + 1;
    
    realold = initSG(minx:maxx, miny:maxy, minz:maxz);
    condi = cond(minx:maxx, miny:maxy, minz:maxz);
    struc = strucmap(minx:maxx, miny:maxy, minz:maxz);
    labelold = labelmap(minx:maxx, miny:maxy, minz:maxz);
    
    %% choose a new patch from TI
    dis_1 = dis3D(TI,realold);
    dis_2 = dis3D(TI,condi);
    dis = (1-w_c)*dis_1 + (w_c*dis_2);
    
    %% candidate that match the conditioning data
    patchcenter = [centerx-minx+1, centery-miny+1, centerz-minz+1];
    canTI = TI(patchcenter(1):TIdim(1)-patchsizex+patchcenter(1),patchcenter(2):TIdim(2)-patchsizey+patchcenter(2),patchcenter(3):TIdim(3)-patchsizez+patchcenter(3));
    candidateind = find(canTI == condi(patchcenter(1),patchcenter(2),patchcenter(3)));
    [candidatex,candidatey,candidatez] = ind2sub(size(dis),candidateind);

    fit = dis(candidateind);
    [fit_v,fit_i] = sort(fit);
    nbm = size(candidatex,1);
    if nbm >= nbreplicates
        fc = randi(nbreplicates);
    else if (nbm<nbreplicates) && (nbm>= 1)
            fc = randi(nbreplicates+1)-1;
        else % nbm == 0
            continue;
        end
    end
    TIpickx = candidatex(fit_i(fc));
    TIpicky = candidatey(fit_i(fc));
    TIpickz = candidatez(fit_i(fc));

    realnew = TI(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1,TIpickz:TIpickz+patchsizez-1);
    labelnew = indexmap(TIpickx:TIpickx+patchsizex-1,TIpicky:TIpicky+patchsizey-1,TIpickz:TIpickz+patchsizez-1);
    
    %% delec boundary
    sigbi = sigb(minx:maxx, miny:maxy, minz:maxz);
    delb = find(sigbi == 1);
    
    %% determine terminal
    errori = abs(realold - realnew);
    c_I = find(isfinite(condi));
    nbci = size(c_I,1);
    
    if nbci > 1
        errold = abs(realold(c_I)-condi(c_I));
        errnew = abs(realnew(c_I) - condi(c_I));
        cstayline = find(errold == 0);
        cchangeline = find(errnew == 0);
        
        Cstay = c_I(cstayline); 
        Cchange = c_I(cchangeline);
        
        T2_c = setdiff(Cchange,Cstay);
        T1_c = setdiff(Cstay, Cchange);
    else
        T2_c = c_I;
        T1_c = [];
    end % devide the conditioning data to two group
    delb = setdiff(delb,T1_c);
    nbT1 = size(T1_c, 1);
    T1ex = [];
    if nbT1*size(T1_c, 2) > 0
        for nt1 = 1:nbT1
            rangmin = errori(T1_c(nt1)) * (1-p);
            rangmax = errori(T1_c(nt1)) * (1+p);
            sigerr = (errori>=rangmin & errori<=rangmax);sigerr=double(sigerr);
            fsig = bwlabeln(sigerr,6);
            T1cn = find(fsig == fsig(T1_c(nt1)));
            T1ex = unique([T1ex;T1cn]);
        end        
    end
    T1ex = setdiff(T1ex,T2_c);
    
%    T1_b = Boundary(patchsizex, patchsizey,patchsizez);
    checkb = ones(patchsizex, patchsizey, patchsizez); 
    checkb(2:patchsizex-1,2:patchsizey-1,2:patchsizez-1)=0;
    T1_b = find(checkb==1);
    T1_b = setdiff(T1_b, T2_c);
    T1_b = setdiff(T1_b, delb);
        
    T1 = union(T1ex, T1_b);
    
    nbT2 = size(T2_c, 1);
    T2ex = [];
    for nt = 1:nbT2
        rangmin = errori(T2_c(nt)) * (1-p);
        rangmax = errori(T2_c(nt)) * (1+p);
        sigerr = (errori>=rangmin & errori<=rangmax);
        fsig = bwlabeln(sigerr,6);
        T2cn = find(fsig == fsig(T2_c(nt)));
        T2ex = unique([T2ex;T2cn]);
    end
    T2 = setdiff(T2_c,T1);
    T2 = setdiff(T2ex,T1);
    
    %% make a new cut
    [oldseamnode,newseamnode,label,preseam] = AddseamCut3D(errori,seamnode,TI, realnew,labelold, labelnew,minx,miny,minz,T1,T2);
    seamnode = newseamnode;
    realprop = realnew;
    realprop(label==0)=realold(label==0);
    initSG(minx:maxx, miny:maxy, minz:maxz)=realprop;
    
    labelprop = labelnew;
    labelprop(label==0) = labelold(label==0);
    labelmap(minx:maxx, miny:maxy, minz:maxz)=labelprop;
    
    ns = ns+1;
    struc(label == 1) = ns;
    strucmap(minx:maxx, miny:maxy, minz:maxz) = struc;
    
    condierr = abs(initSG - cond);
    condierr(isnan(condierr)) = 0;
    a = find(condierr ~= 0);
    nbch = size(a,1);
       
end

initnamei = ['con_cond',filename,num2str(rr),'.SGEMS'];
WriteGrid(initSG,initnamei,'con_cond');
WriteGrid(strucmap,['con_cond-struct-',num2str(rr),'.SGEMS'],'cond');

figure(1); clf;
plot([1:size(meancoerr,1)],meancoerr);
title('mean condi-err');
saveas(gcf,['condierr-',num2str(rr),'.fig']);
condSG = initSG;

%% update costmap
costmap = zeros(SGdim(1:3));
costmap(seamnode(:,5)) = 0.5 * seamnode(:,4);
costmap(seamnode(:,6)) = 0.5 * seamnode(:,4);

end