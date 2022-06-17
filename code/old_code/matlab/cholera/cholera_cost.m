%-------------------------------------------
% cholera_cost.m
% 
% This function encodes the objective functional
% in our cholera model. We implement a left-endpoint
% approximation for the integral.
% 
% Author: George H Lytle
% Date: June 2021
%-------------------------------------------
function [J1b,J1c,J1e,J2b,J2c,J2e] = cholera_cost(oc_params, params, x, v1, v2, tfinal, tstep)
%Translate the parameter vectors
b1 = oc_params(1);          b2 = oc_params(2);
C1 = oc_params(3);          C2 = oc_params(4);
epsilon1 = oc_params(5);    epsilon2 = oc_params(6);
beta_I1 = params(3);        beta_I2 = params(4);
beta_W1 = params(5);        beta_W2 = params(6);

%Translate State and control vectors
S1 = x(:, 1);   I1 = x(:, 2);   W1 = x(:, 4); 
S2 = x(:, 5);   I2 = x(:, 6);   W2 = x(:, 8); 

%Compute total number of subintervals for the simulation
N = tfinal/tstep;

%Compute the objective functional
integral1b = b1.*(beta_I1*S1.*I1 + beta_W1*S1.*W1) ;
integral1c = C1*v1.*S1;
integral1e = epsilon1*v1.^2;
integral2b = b2.*(beta_I2*S2.*I2 + beta_W2*S2.*W2);
integral2c = C2*v2.*S2;
integral2e = epsilon2*v2.^2;
J1b = sum(integral1b(1:N))*tstep;
J1c = sum(integral1c(1:N))*tstep;
J1e = sum(integral1e(1:N))*tstep;
J2b = sum(integral2b(1:N))*tstep;
J2c = sum(integral2c(1:N))*tstep;
J2e = sum(integral2e(1:N))*tstep;

%J = J1+J2;
end

