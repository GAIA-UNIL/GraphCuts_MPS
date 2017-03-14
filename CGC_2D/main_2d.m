% main function
clear; home
clc
addpath 'maxflow';

% This code is for unconditional and conditional patch-based geostatistical
% simulation using graph cuts

% List of Input:
% 1. TI: the source training image
% 2. patchsize=[patchsizex,patchsizey]
% 3. overlapsize=[overlapsizex,overlapsizey]
% 4. nbcandidates: number of candidates
% 5. SGsize=[SGsizex,SGsizey]
% 6. c: conditioning data set
% 7. nc: number of conditioning data
% 8. w: weight for conditioning
% 9. condtype: 0 for unconditional simulation,1 for conditional
% 10. variabletype: 0 for categorical, 1 for continuous variables
% 11. showfig: 0 show figures, 1 don't show the figures
% 12. w_v: Each variable's weight (for multivariable)
% 13. 


%% Define
global nbreplicates patchsize overlap
global maxloop stop
global showfig
global w_v variabletype nbvar
global TIdim SGdim indexmap
global patchrange r w
global th p

%% Initialization
Tiname = 'channel2002.SGEMS';
condname = 'channel_100.SGEMS';
SGname = 'SG';
varname = 'test';

condtype = 1; % 1 for conditioning and 0 for unconditioning
variabletype = [1];
showfig = 1; % 1: show intermediate figures; 0: don't show pictures

nbsim = 1;
nbreplicates = 10;
patchsize=[60,60];
overlap = [15,15];
maxloop = 10;
stop = 0.01;
r = 1.3;
p =0.5;
w_v = [1];
nbvar = 1;
w = 0.8;
th = 0;
% w1 = 0.1; % weight for conditioning data
range = [30 80;30 80; 1 1]; % patch range for x,y,and z for random path

if overlap(1)>patchsize(1) || overlap(2)>patchsize(2)
    disp('wrong input for the first loop');
end
if nbsim<1
    nbsim = 1;
end

% if condtype == 1
%     Condname = 'paola_samp50.SGEMS';
%     [c,nbc,namevar]=LoadPts(Condname);
% end

patchrange = [30,40;30,40];

%% input the training image
TI = LoadGrid(Tiname);
TIdim = size(TI);
SGdim = TIdim;
SG = nan(SGdim);

N_ti = TIdim(1)*TIdim(2);
indexmap = [1:N_ti];
indexmap = reshape(indexmap,TIdim(1),TIdim(2));

%% the main implementation part
if condtype == 0
    GC_MPS_2D(TI,nbsim,SGname);
else
    [conditioning, nbcond, namevar] = LoadPts(condname);
    Cond = nan(SGdim);
    indcond = sub2ind(SGdim(1:2),conditioning(:,1),conditioning(:,2));
    Nsimnode = SGdim(1)*SGdim(2);
    for i = 0:nbvar-1
        Cond(indcond + (i*Nsimnode)) = conditioning(:,4+i);
    end
    CGC_MPS_2D(TI,Cond,nbsim,SGname);
end


%% analysis step
%variocompare_2D(TI,initSG,SGi,SGn,count,validloop);

