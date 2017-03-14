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
%       
%          1
%         ---    cut
%        4 2 5
%          3  
%Template is turned around its axis      
%
%height=5;
%width=5;

N = height*width;
E = [];
% 
%          1
%         ---    cut
%        4 2 5
%          3  
%
is = [1:N]'; is([height:height:N])=[];
js = is+1;
ks = js+1;      ks([height-1:height-1:N-width])=js([height-1:height-1:N-width]);
ls = js-height; ls(1:height-1)=js(1:height-1);
ms = js+height; ms((height-1)*(width-1):N-width)=js((height-1)*(width-1):N-width);
E = [E;is js ks ls ms];

%          3
%        5 2 4
%         ---   cut
%          1  

ks = is-1;       ks([1:height-1:N-width])=is([1:height-1:N-width]);
E = [E;js is ks ms-1 ls-1];
%   
%          5
%        1|2 3
%          4  
is = [1:N-height]';
js = is+height;
ks = js+height; ks(N-2*height+1:N-height)=js(N-2*height+1:N-height);
ls = js+1;      ls(height:height:N-height)=js(height:height:N-height);
ms = js-1;      ms(1:height:N-height)=js(1:height:N-height);
E = [E;is js ks ls ms];
%   
%          4
%        3 2| 1
%          5  
ks = is-height; ks(1:height)=is(1:height);
ls = is-1;      ls(1:height:N-height)=is(1:height:N-height);
ms = is+1;      ms(height:height:N-height)=is(height:height:N-height);
E = [E;js is ks ls ms];

  end