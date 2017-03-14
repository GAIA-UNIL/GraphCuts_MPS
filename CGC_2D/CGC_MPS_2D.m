function CGC_MPS_2D(TI,Cond,nbsim,SGname)

for rr = 1:nbsim
    %% the first step
    [initSG,seamnode,costmap,labelmap,errorco,unsimcond] = init_2D_con(rr,TI,Cond,SGname);
    
    %% the second step: loop for un-matched conditioning data
    [condSG,seamnode,costmap,labelmap,errorcoi] = cond_2D_con(rr,TI,initSG,Cond,SGname,seamnode,errorco,unsimcond,labelmap,costmap);
    
    %% loop for seam nodes
    [modiSG,meancost,seamnode] = modi_2D_con_v2(rr,condSG,TI,SGname,seamnode,Cond,costmap,errorco,labelmap);
end

end