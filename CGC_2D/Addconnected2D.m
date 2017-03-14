function [Connected, Nedge] = Addconnected2D(C)

% ADDCONNECTED creates edges where a node is connected to its
% neighbors that has a different finite values
% Matric C has three possible values nan, 0, and 1
% Connected(:,1) --> Connected(:,2)-->Connected(:,3) the added node index, 
% index of the nodes of value 0,  index of the node of value 1
% Nadd is the number of the added nodes

[height, width] = size(C);
if height<2 && width<2
    disp('No edges are found because the size of the matrix is too small!');
end

N = height * width;
%userview = memory;
%userview
AC = nan(N*4,3);

nb = 0;
no = 0;
Nedge = 0;

% Search down
for is = 1 : N-1
    flag = nan;
    if rem(is, height) ~= 0
        flag = C(is) - C(is + 1);
    end
    if flag == -1
        Nedge = Nedge + 1;
        nodeindex = N + Nedge;
        nb = nb + 1;
        AC(nb,:) = [nodeindex, is, is+1];
    else if flag == 1
             Nedge = Nedge + 1;
             nodeindex = N + Nedge;
             nb = nb + 1;
             AC(nb, :) = [nodeindex, is+1, is];
        end
    end
end

% Search left
for js = height+1:N
    flag = C(js) - C(js - height);
    if flag == -1
        Nedge = Nedge + 1;
        no = no + 1;
        nodeindex = N + Nedge;
        AC(nb+no,:) = [nodeindex, js, js-height];
    else if flag == 1
            Nedge = Nedge + 1;
            no = no + 1;
            nodeindex = N + Nedge;
            AC(nb+no,:) = [nodeindex,js-height,js];
        end
    end
end

Connected = AC(1:nb+no,:);
clear AC;

end