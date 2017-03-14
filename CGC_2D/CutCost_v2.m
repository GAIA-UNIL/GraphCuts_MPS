function cutcost = CutCost_v2(v1,v2)
% This function calculate the difference between two vector v1 and v2

global w_v;
global variabletype;

n = size(w_v,2);
m = size(v1,1);
cutcost = zeros(m,1);
for i = 1:n
   if variabletype(i) == 0
       partcost = v1(:,i) - v2(:,i);
       partcost(partcost ~= 0) = 1;
       cutcost = cutcost + (partcost*w_v(i));
   else
       partcost = abs(v1(:,i) - v2(:,i))*w_v(i);
       cutcost = cutcost + partcost;
   end
end


end