
global nbreplicates patchsize overlap patchrange nbsim
global maxloop stop
global showfig
global w_v variabletype nbvar condtype
global TIdim SGdim
global indexmap TI
global r w th p

condtype = 0; % 1 for conditioning and 0 for unconditioning
variabletype = 0; % 0 for categorical variables
showfig = 0; % 1: show intermediate figures; 0: don't show pictures
nbreplicates = 10;
patchsize=[50,50];  overlap = [15,15];
maxloop = 100;
stop = 0.1;
r = 1.3;    p =0.5;
w_v = 1;  nbvar = 1;  w = 0.8;
th = 0;
nbsim = 100;
patchrange = [20 40;20 40; 1 1]; % patch range for x,y,and z for random path

if overlap(1)>patchsize(1) || overlap(2)>patchsize(2)
    disp('wrong input for the first loop');
end
if nbsim<1
    nbsim = 1;
end

Tiname = 'channel2002.SGEMS';
SGname = 'SG';
varname = 'test';

TI = LoadGrid(Tiname);

TIdim = size(TI);
SGdim = TIdim;

N_ti = TIdim(1)*TIdim(2);
indexmap = [1:N_ti];
indexmap = reshape(indexmap,TIdim(1),TIdim(2));
