function [ y_or_n ] = if_Z_inset(Z,rho)

SZ=cumsum(Z.^2);
rho_v=(1:1:length(Z))*rho;
tv1=(SZ-rho_v);
tv2=(tv1<=0);
 y_or_n=(sum(tv2(:))==length(Z));

end

