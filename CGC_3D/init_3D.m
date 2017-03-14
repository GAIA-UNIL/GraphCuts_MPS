function [initSG,seamnode,costmap,labelmap,strucmap] = init_3D(rr,TI,filename)
% This is the function to generate the initial realization
global nbreplicates patchsize overlap
global indexmap SGdim TIdim

initSG = nan(SGdim);    costmap = zeros(SGdim); labelmap = nan(SGdim);
offset = patchsize - overlap;
strucmap = zeros(SGdim); ns = 0;

%sigb = ones(SGdim); sigb(2:SGdim(1)-1,2:SGdim(2)-1,2:SGdim(3)-1)=0;
sigb = zeros(SGdim); 
sigb(SGdim(1),:,:) = 1; sigb(:,SGdim(2),:) = 1;
seamnode = nan(0,10);

for kk = 1:offset(3):SGdim(3)-overlap(3)
    for jj = 1:offset(2):SGdim(2)-overlap(2)
        for ii = 1:offset(1):SGdim(1)-overlap(1)
            maxx = min(ii+patchsize(1)-1,SGdim(1));
            maxy = min(jj+patchsize(2)-1,SGdim(2));
            maxz = min(kk+patchsize(3)-1,SGdim(3));
            
            simsizex = maxx - ii + 1;
            simsizey = maxy - jj + 1;
            simsizez = maxz - kk + 1;
            
            if ii*jj*kk == 1
                TIpickx = randi(TIdim(1)-simsizex+1);
                TIpicky = randi(TIdim(2)-simsizey+1);
                TIpickz = randi(TIdim(3)-simsizez+1);
                initSG(ii:maxx,jj:maxy,kk:maxz)=TI(TIpickx:TIpickx+simsizex-1,TIpicky:TIpicky+simsizey-1,TIpickz:TIpickz+simsizez-1);
                labelmap(ii:maxx,jj:maxy,kk:maxz)=indexmap(TIpickx:TIpickx+simsizex-1,TIpicky:TIpicky+simsizey-1,TIpickz:TIpickz+simsizez-1);
            else
                ns = ns + 1;
                simold = initSG(ii:maxx,jj:maxy,kk:maxz);
                sigbi = sigb(ii:maxx,jj:maxy,kk:maxz);
                labelold = labelmap(ii:maxx,jj:maxy,kk:maxz);
                struc = strucmap(ii:maxx,jj:maxy,kk:maxz);

                %% find overlap
                s = find(isfinite(simold));
                [sx,sy,sz] = ind2sub([simsizex,simsizey,simsizez],s);
                patchx=max(sx(:)); patchy=max(sy(:)); patchz=max(sz(:));
                realold = simold(1:patchx,1:patchy,1:patchz);
                sigbi = sigbi(1:patchx,1:patchy,1:patchz);

                %% determine termianl
                sigter = zeros(patchx,patchy,patchz);
                if kk==1 && jj==1 && ii~= 1
                    sigter(1,:,:) = 1; % left
                    sigter(patchx,:,:) = 2;
                else if kk==1 && jj~=1 && ii==1
                        sigter(:,1,:) = 1; % bottom
                        sigter(patchx,:,:) = 1; % right
                        sigter(:,patchy,:) = 2; % top
                    else if kk==1 && jj~=1 && ii~=1
                            sigter(1,:,:) = 1; % left
                            sigter(:,1,:) = 1; % bottom
                            sigter(patchx,:,:) = 1; % right
                            sigter(overlap(1):patchx,overlap(2):patchy,:) = 2;
                        else if kk~=1 && jj==1 && ii==1
                                sigter(:,:,1) = 1; % back
                                sigter(:,patchy,:) = 1; % top
                                sigter(patchx,:,:) = 1; % right
                                sigter(:,:,patchz) = 2; % front
                            else if kk ~= 1 && jj==1 && ii~=1
                                    sigter(1,:,:) = 1; sigter(patchx,:,:)=1; % left and right
                                    sigter(:,patchy,1:overlap(3)-1) = 1; % top
                                    sigter(:,:,1) = 1; % back
                                    sigter(overlap(1):patchx,:,overlap(3):patchz) = 2;
                                else if kk~=1 && jj~=1 && ii==1
                                        sigter(:,1,:) = 1; sigter(:,patchy,:)=1; % bottom and top
                                        sigter(:,:,1) = 1; % back
                                        sigter(patchx,:,:) = 1; % right
                                        sigter(:,overlap(2):patchy,overlap(3):patchz) = 2;
                                    else
                                        sigter(1,:,:) = 1; sigter(patchx,:,:) = 1; % left and right
                                        sigter(:,:,1) = 1; % back
                                        sigter(:,patchy,1:overlap(3)-1) = 1; % top
                                        sigter(overlap(1):patchx,overlap(2):patchy,overlap(3):patchz) = 2;
                                    end
                                end
                            end
                        end
                    end
                end  % all cases                                          
                                            
                %% there are three kind of boundaries
                % 1: boundary of simulation grid. 2: outer boundary of the
                % overlap. 3: inter boundary of the overlap. 1 is free, 2 is
                % the source and 3 is the sink

                %% boundray 1:
                sb = find(sigbi == 1);
                
                T1=find(sigter == 1); 
                T1 = setdiff(T1,sb);
                T2=find(sigter==2); 
                
                %% search the training image to find a new patch
                ignx=simsizex-patchx; igny=simsizey-patchy; ignz=simsizez-patchz;
                TIs = TI(1:TIdim(1)-ignx,1:TIdim(2)-igny,1:TIdim(3)-ignz);
                dis = dis3D(TIs,realold);

                %% choose one that has lowest distance
                [v,i]=sort(dis(:));
                sel = randi(nbreplicates);
                ind = i(sel);
                [TIpickx,TIpicky,TIpickz]=ind2sub(size(dis),ind);
                simnew = TI(TIpickx:TIpickx+simsizex-1,TIpicky:TIpicky+simsizey-1,TIpickz:TIpickz+simsizez-1);
                realnew = simnew(1:patchx,1:patchy,1:patchz);
                labelnew = indexmap(TIpickx:TIpickx+simsizex-1,TIpicky:TIpicky+simsizey-1,TIpickz:TIpickz+simsizez-1);
                
                labelrealold = labelold(1:patchx,1:patchy,1:patchz);
                labelrealnew = labelnew(1:patchx,1:patchy,1:patchz);

                realold(isnan(realold)) = realnew(isnan(realold));
                errori = abs(realold - realnew);

                %% make a new cut
                [oldseamnode,newseamnode,label,preseam] = AddseamCut3D(errori,seamnode,TI, realnew,labelrealold, labelrealnew,ii,jj,kk,T1,T2);
                seamnode = newseamnode;
                realprop = realnew;
                realprop(label==0)=realold(label==0);
                simprop = simnew;
                simprop(1:patchx,1:patchy,1:patchz) = realprop;
                initSG(ii:maxx,jj:maxy,kk:maxz) = simprop;
                
                labelrealprop = labelrealnew;
                labelrealprop(label==0) = labelrealold(label==0);
                labelprop = labelnew;
                labelprop(1:patchx,1:patchy,1:patchz)=labelrealprop;
                labelmap(ii:maxx,jj:maxy,kk:maxz) = labelprop;
                
                struci = struc(1:patchx,1:patchy,1:patchz);
                struci(label==1) = ns;
                struc(1:patchx,1:patchy,1:patchz) = struci;
                strucmap(ii:maxx,jj:maxy,kk:maxz) = struc;
                
                %% update costmap
                costmap = zeros(SGdim(1:3));
                costmap(seamnode(:,5)) = 0.5 * seamnode(:,4);
                costmap(seamnode(:,6)) = 0.5 * seamnode(:,4);
                
            end

        end % Loop for x direction
    end % Loop for y direction
end % Loop for z direction

initnamei = ['init',filename,num2str(rr),'.SGEMS'];
WriteGrid(initSG,initnamei,'init');
WriteGrid(strucmap,['init-struct-',num2str(rr),'.SGEMS'],'init');

end