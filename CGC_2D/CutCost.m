function cutcost = CutCost(v1,v2)
% This function calculate the difference between two vector v1 and v2

global w_v;
global variabletype;

n = size(w_v,2);
tm = size(v1,1);
m = tm/2;
cutcost = zeros(m,1);
for i = 1:n
   if variabletype(i) == 0
       partcost = abs(v1(:,i)-v2(:,i));
       partcost(partcost ~= 0) =  w_v(i);
       cutcost = (partcost(1:m)+ partcost(m+1:tm))+cutcost;
   else
       partcost = abs(v1(:,i) - v2(:,i));
       partcost(isnan(partcost)) = 0;
       partcost = partcost * w_v(i);
       cutcost = cutcost + partcost(1:m)+partcost(m+1:tm);
   end
end
cutcost = cutcost + eps;


end