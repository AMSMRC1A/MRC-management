%% Ebola: Optimal Control (main)

% Approximates the optimal vaccination strategy (v* = [v1*, v2*]) for the 
% two-patch Ebola model via the Forward-Backward Sweep Method. 

clear all
close all
clc

%% Initialize

% parameters
parameter_settings
M1 = 0.005; % max vaccination rate for patch 1
M2 = 0.005; % max vaccination rate for patch 2

% time interval
T = 750; % final time
dt = 0.1; % time step
tvec = 0:dt:T; % time interval of T/dt evenly spaced points
M = length(tvec); % length of time interval 

% tolerance
delta = 0.01;

% empty vectors
x      = zeros(M,12); % 12 states
lambda = zeros(M,12); % 12 adjoints
v1     = zeros(M,1);  % control for patch 1
v2     = zeros(M,1);  % control for patch 2

%% Find optimal vaccination strategy

test  = -1;
count = 0;

while(test < 0 && count < 100)
%     figure(1)
%     subplot(1,2,1)
%     plot(tvec,v1)
%     subplot(1,2,2)
%     plot(tvec,v2)
    % update states, adjoints, and controls
    oldx = x;
    oldlambda = lambda;
    oldv1 = v1; 
    oldv2 = v2;
    
    % solve state ODEs
    solx = ode45(@(t,x) ebola_optimal_states(t,x,tvec,v1,v2,par),tvec,y0);
    x = deval(solx,tvec)';

    % solve adjoint ODEs
    sollamb = ode45(@(t,lambda) ebola_optimal_adjoints(t,lambda,tvec,x,v1,v2,par),[T 0],zeros(1,12));
    lambda = deval(sollamb,tvec)';
    
    % pull out states and adjoints needed in v*
    S1 = x(:,1);
    E1 = x(:,2);
    I1 = x(:,3);
    S2 = x(:,7);
    E2 = x(:,8);
    I2 = x(:,9);
    lambda1 = lambda(:,1);
    lambda6 = lambda(:,6);
    lambda7 = lambda(:,7);
    lambda12 = lambda(:,12);
    
    % calculate v*    
    temp1 = (-C1*(S1 + E1) + lambda1.*S1 - lambda6.*S1)/(2*epsilon1);
    temp2 = (-C2*(S2 + E2) + lambda7.*S2 - lambda12.*S2)/(2*epsilon2);
    v11 = min(M1,max(0,temp1));
    v21 = min(M2,max(0,temp2));
    
%     figure(2)
%     subplot(1,2,1)
%     plot(tvec,temp1)
%     subplot(1,2,2)
%     plot(tvec,temp2)
    % update control
    v1 = 0.5*(v11 + oldv1);
    v2 = 0.5*(v21 + oldv2);
    
    % tolerance: |oldu - u|/|u| < delta
    delta*sum(abs(v1))-sum(abs(oldv1-v1))
    [delta*sum(abs(v1))-sum(abs(oldv1-v1)) delta*sum(abs(v2))-sum(abs(oldv2-v2)) delta*sum(abs(x))-sum(abs(oldx-x)) delta*sum(abs(lambda))-sum(abs(oldlambda-lambda))]
    %test = min([delta*sum(abs(v1))-sum(abs(oldv1-v1)) delta*sum(abs(v2))-sum(abs(oldv2-v2)) delta*sum(abs(x))-sum(abs(oldx-x)) delta*sum(abs(lambda))-sum(abs(oldlambda-lambda))])
    test = min([delta*sum(abs(v1))-sum(abs(oldv1-v1)) delta*sum(abs(v2))-sum(abs(oldv2-v2)) ])
    % update count
    count = count + 1
    
end

%% Calculate J

% pull out states and adjoints needed in J
S1 = x(:,1);
E1 = x(:,2);
I1 = x(:,3);
D1 = x(:,5);
S2 = x(:,7);
E2 = x(:,8);
I2 = x(:,9);
D2 = x(:,11);

J1 = sum((b1*(betaI1*S1.*I1 + betaD1*S1.*D1) + b2*(betaI2*S2.*I2 + betaD2*S2.*D2)).*dt) % approx. total cost of new cases
J2 = sum((C1*v1.*(S1 + E1) + epsilon1*v1.^2 + C2*v2.*(S2 + E2) + epsilon2*v2.^2).*dt) % approx. total cost of vaccination

%% Plots

figure(1) % patch 1
subplot(4,2,1); plot(tvec,x(:,1),tvec,x(:,7),'linewidth',2)
subplot(4,2,1);   xlabel('Time')
subplot(4,2,1);   ylabel('S')
subplot(4,2,1);   legend('patch 1','patch 2')

subplot(4,2,2); plot(tvec,x(:,2),tvec,x(:,8),'linewidth',2)
subplot(4,2,2);   xlabel('Time')
subplot(4,2,2);   ylabel('E')

subplot(4,2,3); plot(tvec,x(:,3),tvec,x(:,9),'linewidth',2)
subplot(4,2,3);   xlabel('Time')
subplot(4,2,3);   ylabel('I')

subplot(4,2,4); plot(tvec,x(:,4),tvec,x(:,10),'linewidth',2)
subplot(4,2,4);   xlabel('Time')
subplot(4,2,4);   ylabel('H')

subplot(4,2,5); plot(tvec,x(:,5),tvec,x(:,11),'linewidth',2)
subplot(4,2,5);   xlabel('Time')
subplot(4,2,5);   ylabel('D')

subplot(4,2,6); plot(tvec,x(:,6),tvec,x(:,12),'linewidth',2)
subplot(4,2,6);   xlabel('Time')
subplot(4,2,6);   ylabel('R')

subplot(4,2,7); plot(tvec,v1,tvec,v2,'linewidth',2)
subplot(4,2,7);   xlabel('Time')
subplot(4,2,7);   ylabel('v')




% figure(2) % patch 1
% subplot(4,2,1);plot(tvec,x(:,1))
% subplot(4,2,1);xlabel('Time')
% subplot(4,2,1);ylabel('S1')
% subplot(4,2,2);plot(tvec,x(:,2))
% subplot(4,2,2);xlabel('Time')
% subplot(4,2,2);ylabel('E1')
% subplot(4,2,3);plot(tvec,x(:,3))
% subplot(4,2,3);xlabel('Time')
% subplot(4,2,3);ylabel('I1')
% subplot(4,2,4);plot(tvec,x(:,4))
% subplot(4,2,4);xlabel('Time')
% subplot(4,2,4);ylabel('H1')
% subplot(4,2,5);plot(tvec,x(:,5))
% subplot(4,2,5);xlabel('Time')
% subplot(4,2,5);ylabel('D1')
% subplot(4,2,6);plot(tvec,x(:,6))
% subplot(4,2,6);xlabel('Time')
% subplot(4,2,6);ylabel('R1')
% subplot(4,2,7);plot(tvec,v1)
% subplot(4,2,7);xlabel('Time')
% subplot(4,2,7);ylabel('v1')

% figure(3) % patch 2
% subplot(4,2,1);plot(tvec,x(:,7))
% subplot(4,2,1);xlabel('Time')
% subplot(4,2,1);ylabel('S2')
% subplot(4,2,2);plot(tvec,x(:,8))
% subplot(4,2,2);xlabel('Time')
% subplot(4,2,2);ylabel('E2')
% subplot(4,2,3);plot(tvec,x(:,9))
% subplot(4,2,3);xlabel('Time')
% subplot(4,2,3);ylabel('I2')
% subplot(4,2,4);plot(tvec,x(:,10))
% subplot(4,2,4);xlabel('Time')
% subplot(4,2,4);ylabel('H2')
% subplot(4,2,5);plot(tvec,x(:,11))
% subplot(4,2,5);xlabel('Time')
% subplot(4,2,5);ylabel('D2')
% subplot(4,2,6);plot(tvec,x(:,12))
% subplot(4,2,6);xlabel('Time')
% subplot(4,2,6);ylabel('R2')
% subplot(4,2,7);plot(tvec,v2)
% subplot(4,2,7);xlabel('Time')
% subplot(4,2,7);ylabel('v2')

% figure(4) % adjoints for patch 1
% subplot(4,2,1);plot(tvec,lambda(:,1))
% subplot(4,2,1);xlabel('Time')
% subplot(4,2,1);ylabel('S1')
% subplot(4,2,2);plot(tvec,lambda(:,2))
% subplot(4,2,2);xlabel('Time')
% subplot(4,2,2);ylabel('E1')
% subplot(4,2,3);plot(tvec,lambda(:,3))
% subplot(4,2,3);xlabel('Time')
% subplot(4,2,3);ylabel('I1')
% subplot(4,2,4);plot(tvec,lambda(:,4))
% subplot(4,2,4);xlabel('Time')
% subplot(4,2,4);ylabel('H1')
% subplot(4,2,5);plot(tvec,lambda(:,5))
% subplot(4,2,5);xlabel('Time')
% subplot(4,2,5);ylabel('D1')
% subplot(4,2,6);plot(tvec,lambda(:,6))
% subplot(4,2,6);xlabel('Time')
% subplot(4,2,6);ylabel('R1')

% figure(5) % adjoints for patch 2
% subplot(4,2,1);plot(tvec,lambda(:,7))
% subplot(4,2,1);xlabel('Time')
% subplot(4,2,1);ylabel('S2')
% subplot(4,2,2);plot(tvec,lambda(:,8))
% subplot(4,2,2);xlabel('Time')
% subplot(4,2,2);ylabel('E2')
% subplot(4,2,3);plot(tvec,lambda(:,9))
% subplot(4,2,3);xlabel('Time')
% subplot(4,2,3);ylabel('I2')
% subplot(4,2,4);plot(tvec,lambda(:,10))
% subplot(4,2,4);xlabel('Time')
% subplot(4,2,4);ylabel('H2')
% subplot(4,2,5);plot(tvec,lambda(:,11))
% subplot(4,2,5);xlabel('Time')
% subplot(4,2,5);ylabel('D2')
% subplot(4,2,6);plot(tvec,lambda(:,12))
% subplot(4,2,6);xlabel('Time')
% subplot(4,2,6);ylabel('R2')