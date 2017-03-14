function E = edges6connected3D(simoldstary)

[dimx,dimy,dimz] = size(simoldstary);
Nxy = dimx * dimy;
Ehor = []; Ever = [];
for kk = 1:dimz
    Exy = edges4connected(dimx,dimy);
    Eterm = Exy + (kk-1)*Nxy;
    Ehor = [Ehor;Eterm];
    
    if kk < dimz
        is = [1:Nxy]';
        js = is + Nxy;
        Ez = [is,js;js,is];
        Ezt = Ez + (kk-1)*Nxy;
        Ever = [Ever;Ezt];
    end
end
E = [Ehor;Ever];

val = simoldstary(E(:,1)) + simoldstary(E(:,2));
del = find(isnan(val));
E(del,:) = [];
    
end

