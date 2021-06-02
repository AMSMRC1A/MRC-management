function [ dx ] = statelab7( t,x,tvec,u )

b=0.525;
d=0.5;
c=0.0001;
e=0.5;
g=0.1;
a=0.2;

u=pchip(tvec,u,t);
dx=zeros(4,1);
dx(1)=b*x(4)-(d+c*x(3)+u)*x(1);
dx(2)=c*x(1)*x(3)-(e+d)*x(2);
dx(3)=e*x(2)-(g+a+d)*x(3);
dx(4)=(b-d)*x(4)-a*x(3);

end