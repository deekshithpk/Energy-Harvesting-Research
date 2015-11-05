function [ Y] = near_neigh_projection( X,rho )

alphak=1;
chk=0;
Z=X;
while chk==0
    Z=alphak*Z;
    chk=if_Z_inset(Z,rho);
    alphak=alphak-.01;
end
Y=Z;
end

