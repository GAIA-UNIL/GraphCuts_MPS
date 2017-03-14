function [Connected, Nedge] = Addconnected3D(C)

% ADDCONNECTED creates edges where a node is connected to its
% neighbors that has a different finite values
% Matric C has three possible values nan, 0, and 1
% Connected(:,1) --> Connected(:,2)-->Connected(:,3) the added node index, 
% index of the nodes of value 0,  index of the node of value 1
% Nadd is the number of the added nodes

[height,width,depth] = size(C);
if height<2 && width<2 && depth<2
    disp('No edges are found because the size of the matrix is too small!');
end

N = height * width;
Nt = height * width * depth;
%userview = memory;
%userview
Connected = [];
Nedge = 0;

for i = 1:depth
    base = (i-1)*N;
    % Search down
    for is = 1 : N-1
        flag = nan;
        is_i = is+base;
        if rem(is, height) ~= 0
            flag = C(is_i) - C(is_i + 1);
        end
        if flag == -1
            Nedge = Nedge + 1;
            nodeindex = N + Nedge;
            AC = [nodeindex, is_i, is_i+1];
            Connected = [Connected;AC];
        else if flag == 1
                 Nedge = Nedge + 1;
                 nodeindex = N + Nedge;
                 AC= [nodeindex, is_i+1, is_i];
                 Connected = [Connected;AC];
            end
        end
    end

    % Search left
    for js = height+1:N
        js_i = js+base;
        flag = C(js_i) - C(js_i - height);
        if flag == -1
            Nedge = Nedge + 1;
            nodeindex = N + Nedge;
            AC= [nodeindex, js_i, js_i-height];
            Connected = [Connected;AC];
        else if flag == 1
                Nedge = Nedge + 1;
                nodeindex = N + Nedge;
                AC= [nodeindex,js_i-height,js_i];
                Connected = [Connected;AC];
            end
        end
    end
end

% search depth
for ks = 1:Nt-N
    flag = C(ks)-C(ks+N);
    if flag == -1
        Nedge = Nedge + 1;
        nodeindex = N + Nedge;
        AC=[nodeindex,ks,ks+N];
        Connected = [Connected;AC];
    else if flag == 1
            Nedge = Nedge + 1;
            nodeindex = N + Nedge;
            AC=[nodeindex,ks+N,ks];
            Connected = [Connected;AC];
        end
    end
end
Nedge = size(Connected,1);
clear AC;

end