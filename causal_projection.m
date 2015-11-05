function [ Y] = causal_projection( X,rho )

lX=length(X);
Y=X;
sumx=0;
for loop=1:lX
x=X(loop);
x1=sign(x)*sqrt(min(loop*rho-sumx,x^2));
Y(loop)=x1;
sumx=sumx+x1^2;
end

end

