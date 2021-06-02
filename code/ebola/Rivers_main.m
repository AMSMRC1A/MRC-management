%% Main file for Two-patch model of Ebola 

close all
clear

%% start timer
tic

%% Use parameter_settings to call parameters
Rivers_parameters;

%% Create time vector 
t0 = 0;
tf = 1000;
dt = 0.01;
%nt = tf/dt;
%tspan = linspace(t0, tf, nt);
tspan = t0:dt:tf;

%% solve ODE in Ebola_model.m with ode45()
[t,y] = ode45(@Rivers_model, tspan, y0', [], p);

%% plot solution 

figure()
for i = 1:6
    subplot(2,3,i); plot(t,y(:,i));
    hold on;
end
    
    
% Plot Cumulative Cases
figure()
plot(t,y(:,7))
title('Cumulative Incidence Patch 1')

% % rows
% m = 2;
% % columns
% n = 3;
% 
% subplot(m,n, 1)
% hold on
% plot(t,y(:,1))
% plot(t,y(:,1+6))
% title('S1')
% 
% subplot(m,n, 2)
% hold on
% plot(t,y(:,2))
% plot(t,y(:,2+6))
% title('E1')
% 
% subplot(m,n, 3)
% hold on
% plot(t,y(:,3))
% plot(t,y(:,3+6))
% title('I1')
% 
% subplot(m,n, 4)
% hold on
% plot(t,y(:,4))
% plot(t,y(:,4+6))
% title('H1')
% 
% subplot(m,n, 5)
% hold on
% plot(t,y(:,5))
% plot(t,y(:,5+6))
% title('D1')
% 
% subplot(m,n, 6)
% hold on
% plot(t,y(:,6))
% plot(t,y(:,6+6))
% title('R1')
% 
% % Zoom in on exposed cases
% figure(2)
% plot(t,y(:,2))
% title('Exposed Patch 1')
% 
% % Zoom in on infected cases
% figure(3)
% plot(t,y(:,3))
% title('Infecteds Patch 1')




%% calculate R0
%R0 = betaI1*N1*(alpha1/(alpha1+mu1))*(1/(gammaI1+phi1+deltaI1+mu1))+betaD1*N1*(alpha1/(alpha1+mu1))*(deltaI1/(gammaI1+phi1+deltaI1+mu1))*(1/(xi1))

%% stop timer
toc