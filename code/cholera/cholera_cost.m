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
function [J] = cholera_cost(oc_params, params, x, v1, v2, tfinal, tstep)
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
integral1 = b1.*(beta_I1*S1.*I1 + beta_W1*S1.*W1) + C1*v1.*S1 + epsilon1*v1.^2;
integral2 = b2.*(beta_I2*S2.*I2 + beta_W2*S2.*W2) + C2*v2.*S2 + epsilon2*v2.^2;
J1 = sum(integral1(1:N))*tstep;
J2 = sum(integral2(1:N))*tstep;

J = J1+J2;
end

