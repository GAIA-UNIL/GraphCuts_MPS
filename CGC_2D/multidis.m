function dis = multidis(a,b)

global w_v nbvar variabletype

dim = size(a);
dis = zeros(dim(1:2));
for k = 1:nbvar
    if variabletype(k) == 0
        t_ei = abs(a(:,:,k)-b(:,:,k));
        t_ei(isnan(t_ei)) = 0;
        dis = dis + (w_v(k)*t_ei);
    else
        t_ei = abs(a(:,:,k)-b(:,:,k));
        t_ei(t_ei ~= 0) = 1;
        dis = dis + w_v(k)*t_ei;
    end
end

end