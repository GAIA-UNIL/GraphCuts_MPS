function CGC_MPS_3D(TI,condname,filename)
% This function is to generate unconditional 3D multi-geostatistical
% simulation using graph cuts

global nbsim
[conditioning,nbcond,namevar] = LoadPts(condname);

for rr =1:nbsim
    %% the initial step
    [initSG,seamnode,costmap,labelmap,strucmap,condierr] = init_con_3D(rr,TI,conditioning,filename);
    figure(1); ViewGrid(initSG); view(3);
    title('initial co-simulation');
    figure(2); ViewGrid(strucmap); view(3);
    title('initial co-structures');
    %% the conditional step
    if sum(condierr(:)) ~= 0
        [condSG,seamnode,costmap,labelmap,strucmap] = condi_3D(rr,TI,initSG,conditioning, condierr, filename,seamnode,labelmap,costmap,strucmap);
    end
    %% the modification step
    modi_con_3D(rr,TI,condSG,seamnode,costmap,labelmap,strucmap,filename);
    
end
end