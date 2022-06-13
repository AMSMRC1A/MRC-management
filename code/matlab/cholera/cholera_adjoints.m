%-------------------------------------------
% cholera_adjoints.m
% 
% This function encodes the adjoint equations
% arising from the necessary conditions of 
% Pontryagin's principle for optimal control.
% 
% Author: George H Lytle
% Date: June 2021
%-------------------------------------------
function [dL] =  cholera_adjoints(t, L, x, tvec, v1, v2, params, oc_params)
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

%Translate the parameters for optimal control in oc_params
b1 = oc_params(1);      b2 = oc_params(2);
C1 = oc_params(3);      C2 = oc_params(4);
epsilon1 = oc_params(5);    epsilon2 = oc_params(6);
%Next allocate space for dX

dL = zeros(8, 1);
v1 = pchip(tvec, v1, t);
v2 = pchip(tvec, v2, t);
x = interp1(tvec, x, t);


dL(1) = -(b1*(beta_I1*x(2) + beta_W1*x(4)) + C1*v1 + L(1).*(-beta_I1*x(2)-beta_W1*x(4)-v1-m1) + L(2).*(beta_I1.*x(2) + beta_W1.*x(4)) + L(3).*v1+L(5).*m1);
dL(2) = -(b1*beta_I1*x(1) + L(1).*(mu1-beta_I1*x(1)) + L(2).*(beta_I1*x(1) - (gamma1+mu1+delta1+n1)) + L(3)*gamma1 + L(4)*xi1 + L(6)*n1);
dL(3) = -(L(1)*mu1 - L(3).*(mu1+m1) + L(7)*m1);
dL(4) = -(b1*beta_W1*x(1) - L(1).*beta_W1.*x(1) + L(2).*beta_W1.*x(1) - L(4)*(xi1+rho1) + L(8)*rho1);
dL(5) = -(b2*(beta_I2*x(6) + beta_W2*x(8)) + C2*v2 +  L(1).*m2 + L(5).*(-beta_I2*x(6)-beta_W2*x(8)-v2-m2) + L(6).*(beta_I2.*x(6) + beta_W2.*x(8)) + L(7).*v2);
dL(6) = -(b2*beta_I2*x(5) + L(5).*(mu2-beta_I2*x(5)) + L(6).*(beta_I2*x(5) - (gamma2+mu2+delta2+n2)) + L(7)*gamma2 + L(8)*xi2 + L(2)*n2);
dL(7) = -(L(5)*mu2 - L(7).*(mu2+m2) + L(3)*m2);
dL(8) =  -(b2*beta_W2*x(5) - L(5).*beta_W2.*x(5) + L(6).*beta_W2.*x(5) - L(8)*(xi2+rho2));

end



