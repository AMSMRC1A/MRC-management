%-------------------------------------------
% cholera_states.m
% 
% This function encodes the two-patch SIRW model
% into a state differential system
% 
% Author: George H Lytle
% Date: June 2021
%-------------------------------------------
function [dX] =  cholera_states(t, x, tvec, v1, v2, params)
%Below is the key for the state vector x
%  x(1) = S1;       x(5) = S2;
%  x(2) = I1;       x(6) = I2;
%  x(3) = R1;       x(7) = R2;
%  x(4) = W1;       x(8) = W2;

%Below we translate the parameter vector params
mu1 = params(1);
mu2 = params(2);
beta_I1 = params(3);
beta_I2 = params(4);
beta_W1 = params(5);
beta_W2 = params(6);
m1 = params(7);
m2 = params(8);
n1 = params(9);
n2 = params(10);
gamma1 = params(11);
gamma2 = params(12);
delta1 = params(13);
delta2 = params(14);
xi1 = params(15);
xi2 = params(16);
nu1 = params(17);
nu2 = params(18);
rho1 = params(19);
rho2 = params(20);

%Next allocate space for dX

dX = zeros(8, 1);
v1 = interp1(tvec, v1, t);
v2 = interp1(tvec, v2, t);

dX(1) = mu1*(x(2)+x(3))-beta_I1*x(1).*x(2) - beta_W1*x(1).*x(4) -v1.*x(1) - m1*x(1) + m2*x(5);
dX(2) = beta_I1*x(1).*x(2)+beta_W1*x(1).*x(4) - (gamma1+mu1+delta1+n1)*x(2)+n2*x(6);
dX(3) = gamma1*x(2) - (mu1+m1)*x(3)+v1.*x(1)+m2*x(7);
dX(4) = xi1*x(2) - (xi1+rho1)*x(4);
dX(5) = mu2*(x(6)+x(7))-beta_I2*x(5).*x(6) - beta_W2*x(5).*x(8) -v2.*x(5) - m2*x(5) + m1*x(1);
dX(6) = beta_I2*x(5).*x(6)+beta_W2*x(5).*x(8) - (gamma2+mu2+delta2+n2)*x(6)+n1*x(2);
dX(7) = gamma2*x(6) - (mu2+m2)*x(7)+v2.*x(5)+m1*x(3);
dX(8) = xi2*x(6) - (xi2+rho2)*x(8) + rho1*x(4);

end

