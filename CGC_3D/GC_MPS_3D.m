function GC_MPS_3D(TI,filename)
% This function is to generate unconditional 3D multi-geostatistical
% simulation using graph cuts

global nbsim
for rr = 1:nbsim
    %% the first step
    [initSG,seamnode,costmap,labelmap,strucmap] = init_3D(rr,TI,filename);
    figure(1); ViewGrid(initSG); view(3);
    title('initial simulation');
    figure(2); ViewGrid(strucmap); view(3);
    title('initial structures');
    %% the modify step
%    condSG = initSG;
    [modiSG,labelmap,strucmap,ch] = modi_3D(rr,TI,initSG,seamnode,costmap,labelmap,strucmap,filename);
    
    if size(ch,1) > 1
%     figure(4); VireGrid(modiSG); view(3);
%     title('modified simulation');
%     figure(5); ViewGrid(strucmap); view(3);
%     title('modified structures');
    else
        disp('No new patches are accpected');
    end
end

end