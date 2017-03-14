function [condSG,seamnode,costmap,labelmap,errorcoi] = cond_2D_con(rr,TI,initSG,Cond,filename,seamnode,errorco,unsimcond,labelmap,costmap)
% loop for un-matched conditioning data

global patchrange
global SGdim TIdim
global showfig
global w
global indexmap
global nbvar
global w_v
global variabletype
global nbreplicates
global th p

count = 1;
errorcoi = [];
meanerr_cond = [];
sigb = ones(SGdim(1:2));
sigb(2:SGdim(1)-1,2:SGdim(2)-1) = 0;
nbch = size(unsimcond,1);
condSG = initSG;
errorco = 9e9;

while (nbch> 0) && (errorco>th)
    
    %% choose a un-match conditioning data to modify
    un_match = randi(nbch);
    center = unsimcond(un_match,:);
    
    minx = max(1,center(1)-patchrange(1)); 
    miny = max(1,center(2)-patchrange(2)); 
    maxx = min(center(1)+patchrange(1),SGdim(1));
    maxy = min(center(2)+patchrange(2),SGdim(2));
    
    patchsizex = maxx - minx + 1;
    patchsizey = maxy - miny + 1;
    
    realold = condSG(minx:maxx,miny:maxy,:);
    labelold = labelmap(minx:maxx,miny:maxy);
    condi = Cond(minx:maxx,miny:maxy,:);
    sigbi = sigb(minx:maxx,miny:maxy);
    % delect boundary
    delb = find(sigbi == 1);
    
    %% find a best match from the training image
%     uncon_dis = distance2D(TI,realold);
%     con_dis = distance2D(TI,condi);
%     filtermap = ((1-w)* uncon_dis) + (w*con_dis);
    filtermap = dis_con2D_v2(TI,realold,condi,w);
    %% conditioning candidates
    patchcenter = [center(1)-minx+1, center(2)-miny+1];
    canTI = TI(patchcenter(1):TIdim(1)-patchsizex+patchcenter(1),patchcenter(2):TIdim(2)-patchsizey+patchcenter(2));
    [candidatex,candidatey] = find(canTI == condi(patchcenter(1),patchcenter(2)));
    candidateind = sub2ind(size(filtermap),candidatex,candidatey);
    fit = filtermap(candidateind);
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

    realnew = TI(TIpickx:TIpickx+patchsizex-1, TIpicky:TIpicky+patchsizey-1,:);
    labelnew = indexmap(TIpickx:TIpickx+patchsizex-1, TIpicky:TIpicky+patchsizey-1);
    
    %% determine the terminal
    errori = multidis(realnew,realold);
    c_It = find(isfinite(condi));
    [Icx_i,Icy_i,Ind] = ind2sub(size(realnew),c_It);
    Ic_i = sub2ind([patchsizex,patchsizey],Icx_i,Icy_i);
    Ic_i = unique(Ic_i);
    [fIcx,fIcy] = ind2sub(SGdim(1:2),Ic_i);

    nbci = size(Ic_i,1);
    if nbci > 1
        cerrnew = multidis(realnew,condi);
        cerrold = multidis(realold,condi);
        errold = cerrold(Ic_i);
        errnew = cerrnew(Ic_i);
        cstayline = find(errold == 0);
        changeline = find(errnew == 0);
        
        cstay = Ic_i(cstayline);
        cchange = Ic_i(changeline);
        
        T2_c = setdiff(cchange,cstay);
        T1_c = setdiff(cstay,cchange);
    else
        T2_c = Ic_i;
        T1_c = [];
    end
    delb = setdiff(delb,T1_c);
    nbT1 = size(T1_c,1);
    T1ex = [];
    
    if nbT1*size(T1_c,2) > 0
        for nt1 = 1:nbT1
            rangmin = errori(T1_c(nt1)) * (1-p);
            rangmax = errori(T1_c(nt1)) * (1+p);
            sigerr = (errori>=rangmin & errori<=rangmax);
            fsig = bwlabel(sigerr);
            T1cn = find(fsig == fsig(T1_c(nt1)));
            T1ex = unique([T1ex;T1cn]);
        end
    end
    
    T1ex = setdiff(T1ex,T2_c);
    
    TB = ones(patchsizex,patchsizey);
    TB(2:patchsizex-1,2:patchsizey-1) = 0;
    
    T1_b = find(TB==1);
    T1_b = setdiff(T1_b, T2_c);
    T1_b = setdiff(T1_b, delb);
        
    T1 = union(T1ex, T1_b);
    
    nbT2 = size(T2_c, 1);
    T2ex = [];
    if nbT2*size(T2_c,2) > 0
        for nt = 1:nbT2
            rangmin = errori(T2_c(nt)) * (1-p);
            rangmax = errori(T2_c(nt)) * (1+p);
            sigerr = (errori>=rangmin & errori<=rangmax);
            fsig = bwlabel(sigerr);
            T2cn = find(fsig == fsig(T2_c(nt)));
            T2ex = unique([T2ex;T2cn]);
        end
    end
    T2 = setdiff(T2_c,T1);
    T2 = setdiff(T2ex,T1);
    
    [oldseamnode,newseamnode,label,preseam] = AddseamCut2D(errori,seamnode,TI, realnew,labelold, labelnew,minx,miny,T1,T2);

    %% update
    seamnode = newseamnode;
    Np = patchsizex*patchsizey;
    I = find(label==0);
    for k = 1:nbvar
        AddI = I + (k-1)*Np;
        I = [I;AddI];
    end

    realprop = realnew;
    realprop(I) = realold(I);
    condSG(minx:maxx,miny:maxy,:) = realprop;

    labelprop = labelnew;
    labelprop(label == 0) = labelold(label == 0);

    labelmap(minx:maxx,miny:maxy) = labelprop;

    %% update costmap
    costmap = zeros(SGdim(1:2));
    costmap(seamnode(:,5)) = 0.5 * seamnode(:,4);
    costmap(seamnode(:,6)) = 0.5 * seamnode(:,4);

    %% update conditioning data
    Conddif = abs(condSG - Cond);
    Conddif(isnan(Conddif)) = 0;
    Ic = find(isfinite(Cond));
    coerrmap = zeros(SGdim(1:2)); % the map for conditioning error
    errorco = 0;
    for k = 1:nbvar
        if variabletype(k) == 0
            Conddif(Conddif(:,:,k)~=0) = 1;                    
        end
        coerrmap = (Conddif(:,:,k) * w_v(k)) + coerrmap;
    end
    errorco = sum(coerrmap(:));
    [Icx,Icy] = ind2sub(SGdim,Ic);
    unsimcond = [Icx,Icy,ones(size(Icx)),coerrmap(Ic)];
    simedline = find(unsimcond(:,4) <= th);
    unsimcond(simedline,:) = [];

    nbch = size(unsimcond,1);
    errorcoi = [errorcoi;errorco];
 %   disp(nbch);
    if showfig == 1
        Td = zeros(patchsizex,patchsizey);
        Td(T1) = 1;    Td(T2) = 2;
        showfigures(label,Td,preseam,costmap,condSG,rr,count,'_cond');
    end
    count = count + 1;
end    


end