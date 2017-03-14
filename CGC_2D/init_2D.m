
function [initSG, seamnode, costmap,labelmap,patchpick] = init_2D(rr,TI,outputname)

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

initSG = nan(SGdim);
costmap = zeros(SGdim(1),SGdim(2));
offset = patchsize - overlap;
labelmap = nan(SGdim(1),SGdim(2));
patchpick = [];
seamnode = [];
for ii = 1:offset(1):SGdim(1)-overlap(1)
    for jj = 1:offset(2):SGdim(2)-overlap(2)
        %% define the real simulation patch
        maxx = min(ii+patchsize(1)-1, SGdim(1));
        maxy = min(jj+patchsize(2)-1, SGdim(2));
        
        simsizex = maxx-ii+1;
        simsizey = maxy-jj+1;
        simold = initSG(ii:maxx, jj:maxy,:); % old patch
        labelold = labelmap(ii:maxx, jj:maxy);
        
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
        realoldstay = realold;
        
        %% find the best matches from training image
        if flag == 0
            TIpickx = randi(TIdim(1)-simsizex);
            TIpicky = randi(TIdim(2)-simsizey);
            initSG(1:simsizex,1:simsizey,:) = TI(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1,:);
            labelmap(1:simsizex,1:simsizey) = indexmap(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1);
            step = 0;
        else
            T1 = find(Td == 1);
            T2 = find(Td == 2);
            [filtermap,No] = distance2D_v1(TI,simold);
            [v,i] = sort(filtermap(:));
            if nbreplicates > 1
                sel = i(randi(nbreplicates));
            else
                sel = i(1);
            end
            [fx,fy] = size(filtermap);
            [TIpickx,TIpicky] = ind2sub([fx,fy],sel);
            simnew = TI(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1,:);
            realnew = simnew(1:px,1:py,:);
            simprop = simnew;
            errori = zeros(px,py);
            
            labelnew = indexmap(TIpickx:TIpickx+simsizex-1, TIpicky:TIpicky+simsizey-1);
            labelrealnew = labelnew(1:px,1:py);
            labelprop = labelnew;
            
            realold(isnan(realold)) = realnew(isnan(realold));
            for k = 1:nbvar
                if variabletype(k) == 1
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
%            [oldseamnode,newseamnode,label,preseam] = AddseamCut2D_v2(errori,seamnode,TI, realnew,labelrealold, labelrealnew,ii,jj,T1,T2);

%             if flag == 1
%                 E = edges4connected(px,py);
%                 V = errori(E(:,1)) + errori(E(:,2)) + eps;
%                 A = sparse(E(:,1),E(:,2),V,Np,Np,4*Np);
%                 size1 = size(T1,1); size2 = size(T2,1);
%                 T = sparse([T1;T2],[ones(size1,1);ones(size2,1)*2],ones(size1+size2,1)*9e9, Np, 2);
%                 [flow,labels] = maxflow(A,T);
%                 label = reshape(labels,px,py);
% 
%                 Nadd = 0; lineno = []; check = []; preseam = [];
%             else
%                 %% Scan the seamlist to find previous seam nodes and make a new cut
%                 [lineno,flowold,flow,labels,Nadd,check,preseam] = Addseam2D(errori,seamnode,TI,realnew,ii,jj,T1,T2);
%                 
%                 %% check new cuts
%                 label = reshape(labels(1:Np),px,py);
% 
%             end % make a new cut

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
            labelprop(1:px,1:py) = labelrealprop;
            labelmap(ii:maxx, jj:maxy) = labelprop;
            
            %% cudate costmap
            costmap = zeros(SGdim(1:2));
            costmap(seamnode(:,5)) = 0.5 * seamnode(:,4);
            costmap(seamnode(:,6)) = 0.5 * seamnode(:,4);

%             %% connection between cuts
%             [AC,Ne] = Addconnected2D(label);
% 
%             %% new seam nodes
%             newseamnode = newseam(AC,ii,jj,1,labelrealnew,labelrealold, errori,TI);
% 
%             %% update seam nodes
%             [seamnode,costmap] = updateseam2D(seamnode,newseamnode,costmap,Nadd,lineno,check);
            
            if showfig == 1
                showfigures(label,Td,preseam,costmap,initSG,rr,step,'_init');
            end
        end
        patchpicki = [simsizex,simsizey,TIpickx,TIpicky,ii,jj];
        patchpick = [patchpick;patchpicki];
        step = step + 1;

    end % Loop for direction y
end % Loop for direction x

initnamei = ['init',outputname,num2str(rr),'.SGEMS'];
WriteGrid(initSG,initnamei,'init');
end