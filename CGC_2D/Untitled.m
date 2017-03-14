% test graph cut
clear all;

%% parameter
global SGdim
global nbvar
global variabletype
global w_v
SGdim = [8 8];
nbvar = 1;
variabletype = 1;
w_v = 1;

vmin = 0;   vmax = 8;
%%
seamlist = [];

% a = LoadGrid('a.txt');
% b = LoadGrid('b.txt');
% c = LoadGrid('c.txt');
a=[]; b=[];

for i = 1:8
%     a1 = i*zeros(1,32);
%     b1 = i*zeros(32,1);
    a1 = [1:8];
    b1 = [1:8]';
    a = [a;a1];
    b = [b,b1];
end
c = 2*ones(1,8);
c = [c;4*ones(1,8)];
c = [c;6*ones(1,8)];
c = [c;8*ones(1,8)];
c = [c;7*ones(1,8)];
c = [c;5*ones(1,8)];
c = [c;3*ones(1,8)];
c = [c;1*ones(1,8)];
figure(3); clf;
imagesc(c); axis xy equal tight; caxis([vmin vmax]);

TI = [a,b,c];

figure(1); clf;
imagesc(a); 
axis xy equal tight; caxis([vmin vmax]);
figure(2); clf;
imagesc(b); axis xy equal tight; caxis([vmin vmax]);

Nindex = 8*24;
N = 64;
index = [1:Nindex]; index = reshape(index,8,24);
T1 = [1:8]';
T2 = [57:64]';

error1 = abs(a-b) + eps;
E = edges4connected(8,8);
V1 = error1(E(:,1)) + error1(E(:,2)) + eps;
A = sparse(E(:,1),E(:,2),V1,N,N,4*N);
size1 = size(T1,1); size2 = size(T2,1);
T = sparse([T1;T2],[ones(size1,1);ones(size2,1)*2],ones(size1+size2,1)*9e9,N,2);
[flow,labels] = maxflow(A,T);
label = reshape(labels,8,8);

labelold = index(1:8,1:8);
labelnew = index(:,9:16);
labelprop = labelnew;
labelprop(label==0)=labelold(label==0);

r1 = b; 
r1(label == 0) = a(label == 0);
[AC,Ne] = Addconnected2D(label);
newseamnode = newseam(AC,1,1,1,labelnew,labelold,error1,[8,8,1]);
costmap = zeros(8,8);
Nadd=0; lineno=[]; check=[];
[seamlist,costmap] = updateseam2D(seamlist,newseamnode,costmap,Nadd,lineno,check);

figure(4); clf;
subplot(2,2,1);
imagesc(error1); axis xy equal tight; caxis([vmin vmax]);
subplot(2,2,2);
imagesc(label); axis xy equal tight; caxis([vmin vmax]);
subplot(2,2,3);
%imagesc(costmap); axis xy equal tight; caxis([0 vmax]);
scatter(seamlist(:,2),seamlist(:,1),20,seamlist(:,4),'filled');
axis xy equal tight; caxis([0 vmax]);
axis([1 8 1 8]);
subplot(2,2,4);
imagesc(r1); axis xy equal tight; caxis([vmin vmax]);

%% The second step
error2 = abs(r1-c)+eps;
indexnew = index(:,17:24);
indexold = labelprop;
[lineno,flowold,flow,labels2,Nadd,check,preseam] = Addseam2D(error2,seamlist,TI,c,1,1,T1,T2);
label2 = reshape(labels2(1:N),8,8);
r2 = c;
r2(label2==0) = r1(label2==0);
[AC,Ne] = Addconnected2D(label2);
newseamnode = newseam(AC,1,1,1,indexnew,indexold,error2,[8,8,1]);
[seamlist,costmap] = updateseam2D(seamlist,newseamnode,costmap,Nadd,lineno,check);

figure(6); clf;
subplot(2,2,1);
imagesc(error2); axis xy equal tight; caxis([vmin vmax]);
subplot(2,2,2);
imagesc(label2); axis xy equal tight; caxis([vmin vmax]);
subplot(2,2,3);
%imagesc(costmap); axis xy equal tight; caxis([0 vmax]);
scatter(seamlist(:,2),seamlist(:,1),20,seamlist(:,4),'filled'); 
axis xy equal tight; caxis([0 vmax]);
axis([1 8 1 8]);
subplot(2,2,4);
imagesc(r2); axis xy equal tight; caxis([vmin vmax]);

