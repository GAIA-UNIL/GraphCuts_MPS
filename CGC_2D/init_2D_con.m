function [initSG,seamnode,costmap,labelmap,errorco,unsimcond] = init_2D_con(rr,TI,Cond,filename)
% The function to generate the initial realization for conditional graph
% cuts

global nbreplicates
global patchsize
global overlap
global showfig
global indexmap
global SGdim
global TIdim
global w_v
global nbvar
global variabletype
global w
global th

%initSG = nan(SGdim);
initSG = Cond;
costmap = zeros(SGdim(1),SGdim(2));
offset = patchsize - overlap;
labelmap = nan(SGdim(1),SGdim(2));
patchpick = [];
seamnode = [];
Ict = find(isfinite(initSG));
[Icx,Icy,Icz,Ind] = ind2sub(SGdim,Ict);
Ic = sub2ind(SGdim(1:2),Icx,Icy);
Ic = unique(Ic);
[fIcx,fIcy] = ind2sub(SGdim(1:2),Ic);

for ii = 1:offset(1):SGdim(1)-overlap(1)
    for jj = 1:offset(2):SGdim(2)-overlap(2)
        %% define the real simulation patch
        maxx = min(ii+patchsize(1)-1, SGdim(1));
        maxy = min(jj+patchsize(2)-1, SGdim(2));
        
        simsizex = maxx-ii+1;
        simsizey = maxy-jj+1;
        simold = initSG(ii:maxx, jj:maxy,:); % old patch
        labelold = labelmap(ii:maxx, jj:maxy);
        
        sigforcondi = Cond(ii:maxx, jj:maxy,:);
        cip = find(isfinite(sigforcondi));
        [pcx,pcy,pcz,pcv] = ind2sub(size(simold),cip);
        pcin = sub2ind([simsizex,simsizey],pcx,pcy);
        TTc = unique(pcin);
        
        %% find a best match from the training image
%         if ii*jj == 1
%             TIpickx = randi(TIdim(1)-simsizex+1);
%             TIpicky = randi(TIdim(2)-simsizey+1);
%         else
%             uncon_dis = distance2D(TI,simold);
%             con_dis = distance2D(TI,sigforcondi);
%             filtermap = ((1-w)* uncon_dis) + (w*con_dis);
% 
%             [v,i] = sort(filtermap(:));
%             if nbreplicates > 1
%                 seli = randi(nbreplicates);
%                 sel = i(seli);
%             else
%                 sel = i(1);
%             end
%             [fx,fy] = size(filtermap);
%             [TIpickx,TIpicky] = ind2sub([fx,fy],sel);
%         end
        filtermap = dis_con2D_v2(TI,simold,sigforcondi,w);
        [v,i] = sort(filtermap(:));
        if nbreplicates > 1
            seli = randi(nbreplicates);
            sel = i(seli);
        else
            sel = i(1);
        end
        [fx,fy] = size(filtermap);
        [TIpickx,TIpicky] = ind2sub([fx,fy],sel);

        simnew = TI(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1,:);

        if ii == 1
            if jj == 1
                flag = 0;
                px = simsizex;
                py = simsizey;
            else
                flag = 1;
                px = simsizex;
                py = overlap(2);
                Np = px * py;
                Td = zeros(px,py);
                Td(:,1) = 1;
                Td(:,py) = 2;
            end
        else
            if jj == 1
                flag = 2;
                px = overlap(1);
                py = simsizey;
                Np = px*py;
                Td = zeros(px,py);
                Td(1,:) = 1;
                Td(:,py) = 1;
                Td(px,1:py-1) = 2;
            else
                flag = 3;
                px = simsizex;
                py = simsizey;
                Np = px*py;
                Td = zeros(px,py);
                Td(:,1) = 1;    Td(1,:) = 1;
                Td(1:overlap(1)-1,py) = 1;
                Td(overlap(1):px,overlap(2):py) = 2;
%                Td(px,1:overlap(2))= 2;
%                 Td(overlap(1),overlap(2)+1:py) = 2;
%                 Td(overlap(1):px,overlap(2)) = 2;
%                 Td(isnan(simold(:,:,1))) = 2;
            end   
        end % different cases for simulation
        %% The working patch
        realold = simold(1:px, 1:py,:);
        labelrealold = labelold(1:px, 1:py);
 %       realoldstay = realold;
        labelnew = indexmap(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1);
        labelrealnew = labelnew(1:px,1:py);

        %% find the best matches from training image
        if flag == 0
            TIpickx = randi(TIdim(1)-simsizex);
            TIpicky = randi(TIdim(2)-simsizey);
            initSG(1:simsizex,1:simsizey,:) = TI(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1,:);
            labelmap(1:simsizex,1:simsizey) = indexmap(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1);
            step = 0;
        else
            T1_b = find(Td == 1);
            T2_b = find(Td == 2);
            condi = sigforcondi(1:px,1:py);
%            
%            Tc = find(isfinite(condi));
            Tc = []; % here we don't need the 
            T1 = union(T1_b,Tc);
            T2 = setdiff(T2_b,Tc);

            realnew = simnew(1:px,1:py,:);
            simprop = simnew;
            errori = zeros(px,py);

%            realold(isnan(realold)) = realnew(isnan(realold));
            for k = 1:nbvar
                if variabletype(k) == 0
                    t_ei = abs(realnew(:,:,k)-realold(:,:,k));
                    t_ei(isnan(t_ei)) = 0;
                    errori = errori + (w_v(k)*t_ei);
                else
                    t_ei = abs(realnew(:,:,k)-realold(:,:,k));
                    t_ei(t_ei ~= 0) = 1;
                    errori = errori + w_v(k)*t_ei;
                end
            end
            
            %% add seam node and make a new cut
            [oldseamnode,newseamnode,label,preseam] = AddseamCut2D(errori,seamnode,TI, realnew,labelrealold, labelrealnew,ii,jj,T1,T2);

            %% update
            seamnode = newseamnode;
            I = find(label==0);
            for k = 1:nbvar
                AddI = I + (k-1)*Np;
                I = [I;AddI];
            end
            
            realprop = realnew;
            realprop(I) = realold(I);

            simprop(1:px,1:py,:) = realprop;
            initSG(ii:maxx, jj:maxy,:) = simprop;

            labelrealprop = labelrealnew;
            labelrealprop(label == 0) = labelrealold(label == 0);
            
            labelprop = labelnew;
            labelprop(1:px,1:py) = labelrealprop;
            labelmap(ii:maxx, jj:maxy) = labelprop;
            
            %% cudate costmap
            costmap = zeros(SGdim(1:2));
            costmap(seamnode(:,5)) = 0.5 * seamnode(:,4);
            costmap(seamnode(:,6)) = 0.5 * seamnode(:,4);
            
            if showfig == 1
                showfigures(label,Td,preseam,costmap,initSG,rr,step,'_init');
            end
        end
        patchpicki = [simsizex,simsizey,TIpickx,TIpicky,ii,jj];
        patchpick = [patchpick;patchpicki];
        step = step + 1;

    end % Loop for direction y
end % Loop for direction x

 %% update error for conditioning data
    Conddif = abs(initSG - Cond);
    Conddif(isnan(Conddif)) = 0;
    coerrmap = zeros(SGdim(1:2));
    errorco = 0;
    for k = 1:nbvar
        if variabletype(k) == 0
            Conddif(Conddif(:,:,k)~=0) = 1;                    
        end
        coerrmap = (Conddif(:,:,k) * w_v(k)) + coerrmap;
    end
    errorco = sum(coerrmap(:));
    unsimcond = [fIcx,fIcy,ones(size(fIcx)),coerrmap(Ic)];
    simedline = find(unsimcond(:,4) <= th);
    unsimcond(simedline,:) = [];

end