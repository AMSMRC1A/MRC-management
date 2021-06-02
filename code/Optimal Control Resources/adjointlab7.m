function [ dlambda ] = adjointlab7( t,lambda,tvec,x,u )

b=0.525;
d=0.5;
c=0.0001;
e=0.5;
g=0.1;
a=0.2;
A=0.1;

x=interp1(tvec,x,t);
u=pchip(tvec,u,t);
dlambda=zeros(4,1);

dlambda(1) = lambda(1)*(d+c*x(3)+u) - c*lambda(2)*x(3);
dlambda(2) = lambda(2)*(e+d) - lambda(3)*e;
dlambda(3) = -A + (lambda(1)-lambda(2))*c*x(1) + lambda(3)*(g+a+d) + lambda(4)*a;
dlambda(4) = -lambda(1)*b - lambda(4)*(b-d);

end

