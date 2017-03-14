function filtermap = distance2D(TI,simold)
% distance between two patches
global w_v;
global nbvar;
global variabletype
global TIdim

row = size(simold,1);
col = size(simold,2);
%filter = ones(row,col);
filtermap = zeros(TIdim(1)-row+1,TIdim(2)-col+1);

for i = 1:nbvar
    A = TI(:,:,i);
    B = simold(:,:,i);
%    filter(isnan(B)) = 0;
%    filtermap(isnan(B)) = 0;

    if variabletype(i) == 0
        % categotical data
        term = zeros(TIdim(1)-row+1,TIdim(2)-col+1);
        exfa = unique(B);
        nanI = find(isnan(exfa));
        exfa(nanI) = [];
        
        nbfa = size(exfa,1);
        for j = 1:nbfa
            Bfi = (B == exfa(j));
            Bfi = single(Bfi);
            Afi = (A == exfa(j));
            Afi = single(Afi);
            filter = ones(size(B));
            filter(B ~= exfa(j)) = 0;
            a2 = filter2(filter,Afi.^2,'valid');
%            b2 = filter2(filter,Bfi.^2,'valid');
            b2 = sum(sum(Bfi.^2));
            ab = filter2(Bfi,Afi,'valid')*2;
            term_k = (a2-ab)+ b2;
            term = term + term_k;
        end
        termd = w_v(i) * sqrt(term)/sum(filter(:));
    else
        B(isnan(B)) = 0;
        filter = ones(size(B));
        filter(isnan(B)) = 0;
        a2 = filter2(filter,A.^2,'valid');
%        b1 = filter2(filter,B.^2,'valid');
        b2 = sum(sum(B.^2));
        ab = filter2(B,A,'valid')*2;
        term = (a2-ab) + b2;
        termd = w_v(i) * sqrt(term)/sum(filter(:));
    end
    filtermap = filtermap + termd;
end
end