% the main function
clear; home
clc
addpath 'maxflow'

% This code is for 3D patch-based geostatistical simulation using graph
% cuts

% List of Input:
% 1. TI: the source training image
% 2. patchsize=[patchsizex,patchsizey]
% 3. overlapsize=[overlapsizex,overlapsizey]
% 4. nbreplicates: number of candidates
% 5. SGsize=[SGsizex,SGsizey]
% 6. c: conditioning data set
% 7. nc: number of conditioning data
% 8. w: weight for conditioning
% 9. condtype: 0 for unconditional simulation,1 for conditional
% 10. variabletype: 0 for categorical, 1 for continuous variables
% 11. showfig: 0 show figures, 1 don't show the figures
% 12. w_v: Each variable's weight (for multivariable)

%% Define
global nbreplicates patchsize overlap halfpr 
global TIdim SGdim indexmap nbsim
global w_c condtype variabletype
global showfig stop maxloop

%% Initialization
%TIname = 'Maules_Creek_3D.SGEMS';
%TIname = '3D_tipart.SGEMS';
TIname = 'new3Dti.sgems';
condname = '3Dcon_10.SGEMS';
SGname = 'SG';
varname = 'test';

condtype = 0;   variabletype = 0;
showfig = 0;
nbsim = 1; nbreplicates=10;
patchsize = [20,20,20]; overlap = [8,8,8];
halfpr = [10,10,8];

maxloop = 2; stop = 0.15;

%% Input the training image
TI = LoadGrid(TIname); 
% TI = TI(1:100,1:100,1:80);
TIdim = size(TI);   SGdim = [50,50,40];

N_ti = TIdim(1)*TIdim(2)*TIdim(3);
indexmap =  [1:N_ti];
indexmap = reshape(indexmap,TIdim(1),TIdim(2),TIdim(3));

%% the main implementation part
if condtype == 0
    GC_MPS_3D(TI,SGname);
else
    w_c = 0.8;
    CGC_MPS_3D(TI,condname,SGname);
end

