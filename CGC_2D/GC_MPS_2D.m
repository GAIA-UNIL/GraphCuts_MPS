function GC_MPS_2D(TI,nbsim,filename)
% This function is to generate unconditional realizations using graph cuts
for rr = 1:nbsim
    %% the first step
    [initSG,seamnode,costmap,labelmap] = init_2D(rr,TI,filename);
    
    %% the second step -- modification
    [modiSG,meancost,seamnode] = modi_2D(rr,initSG,TI,filename,seamnode,costmap,labelmap);
    
    %% statistics
%     figure(5); clf;
%     nbloop = size(meancost,1);
%     plot([1:nbloop],meancost);
end

end