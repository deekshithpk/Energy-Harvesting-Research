function [ nx ] = Uniform_nsphere( radis,n )

nx=normrnd(0,radis,1,n);
nxs=nx.^2;
nx=nx/sqrt(sum(nxs(:)));
U=radis*(unifrnd(0,1))^(1/3);
nx=nx*U;


end

