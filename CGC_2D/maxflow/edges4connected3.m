  function E = edges4connected2(height,width)

% EDGE4CONNECTED Creates edges where each node
%   is connected to its four adjacent neighbors on a 
%   height x width grid.
%   E - a vector in which each row i represents an edge
%   E(i,1) --> E(i,2). The edges are listed is in the following 
%   neighbor order: down,up,right,left, where nodes 
%   indices are taken column-major.
%
%   (c) 2008 Michael Rubinstein, WDI R&D and IDC
%   $Revision$
%   $Date$
%   
% number indicates position of element in the row
%        D C B
%          A
%         ---    cut
%          1
%        2 3 4
%1. as shown
%2. turned 180?
%3. turned 90? counterclockwise
%4. turned 90? clockwise
%Template is turned around its axis 
% colums: 1 2 3 4 A B C D 
%
%height=5;
%width=5;

N = height*width;
E = [];
 
as   = [1:N]';       as([height:height:N])=[];
cs   = as-1;         cs(1:(height-1):(N-width))=as(1:(height-1):(N-width));
bs   = cs+height;    bs(((height-1)*(width-1)+1):(N-width))=cs(((height-1)*(width-1)+1):(N-width));
ds   = cs-height;    ds(1:(height-1))=cs(1:(height-1));

is   = as+1;         
iiis = is+1;         iiis((height-1):(height-1):(N-width))=is((height-1):(height-1):(N-width));    
iis  = iiis-height;  iis(1:(height-1))=iiis(1:(height-1));
vis  = iiis+height;  vis(((height-1)*(width-1)+1):(N-width))=iiis(((height-1)*(width-1)+1):(N-width));

E = [E;is iis iiis vis as bs cs ds];
E = [E;as bs cs ds is iis iiis vis];


as = [1:N-height]';
cs = as-height; cs(1:height)=as(1:height);      
bs = cs-1;      bs(1:height:N-height)=cs(1:height:N-height);     
ds = cs+1;      ds(height:height:N-height)=cs(height:height:N-height);    

is = as+height;
iiis = is+height; iiis((N-2*height+1):(N-height))=is((N-2*height+1):(N-height));   
iis = iiis+1;     iis(height:height:N-height)=iiis(height:height:N-height);      
vis = iiis-1;     vis(1:height:N-height)=iiis(1:height:N-height);

E = [E; is iis iiis vis as bs cs ds];
E = [E; as bs cs ds is iis iiis vis];


  %end