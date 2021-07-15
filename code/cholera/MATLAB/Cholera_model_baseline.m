%%%%%%%%%%% ODE solution for Cholera SIR two-patch model;%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  May 31st, 2021 %%%%%%%%%%%


% Parameter values 
% Patch 1
mu1 = 1E-4; 
beta1I=2.64E-5;
beta1W= 1.21E-4;
m1=.1; % S,R movement
delta1=5E-4;
gamma1=.25;
n1=0.001; % I movement
xi1=7.56E-3;
nu1=xi1; % decay of pathogen
rho1=0.1; % pathogen transport rate
% Patch 2
mu2=mu1; 
beta2I=beta1I;
beta2W=beta1W;
m2=m1;
delta2=delta1;
gamma2=.25;
n2=n1;
xi2=xi1;
nu2=xi1;
rho2=.01;


%Reproduction number (ODE with constant transmission rates)
% F=[betAbar*SstA0/NA,0; 0,betBbar*SstB0/NB];
% V=[muA+gamaA+phiIAB, -phiIBA;-phiIAB, muB+gamaB+phiIBA];
% M=F*V^(-1);
% e=eig(M)
% % R0A=;%%expression is too long
% % R0B=;  
% R0=max(e) 




% Timespan in days
tspan=0:1:365;

% Initial conditions
%Patch 1
N10=1e6;
x10=N10;%S1
x20=0;%I1
x30=0;%R1
x40=0;%W1
%Patch 2
N20=N10;
x50=N20;%S2
x60=1;%I2
x70=0;%R2
x80=0;%W2
ICs = [x10, x20, x30, x40, x50, x60, x70, x80];


% Get solutionss
options = odeset('RelTol',10^-8);
[t y]=ode45(@cholera, tspan, ICs, options);
 x1=y(:,1);
 x2=y(:,2);
 x3=y(:,3);
 x4=y(:,4);
 x5=y(:,5);
 x6=y(:,6);
 x7=y(:,7);
 x8=y(:,8);


figure(1)
subplot(2,2,1)
plot( t, [x1,x5], 'linewidth',2)
legend('S1','S2')
subplot(2,2,2)
plot( t, [x2,x6], 'linewidth',2)
legend('I1','I2')
subplot(2,2,3)
plot( t, [x3,x7], 'linewidth',2)
legend('R1','R2')
subplot(2,2,4)
plot( t, [x4,x8], 'linewidth',2)
legend('W1','W2')

% figure(1)
% plot(t, log(x1),'b', 'linewidth',2)
% legend('IA')
% xlim([20,25]);
% figure(1)
% plot(t, x1,'b', 'linewidth',2)
% legend('IA')

% figure(2)
% plot(t, x2,'r', 'linewidth',2)
% legend('IB')
% figure(3)
% plot(t, x5,'b', 'linewidth',2)
% legend('RA')
% figure(4)
% plot(t, x6,'b--', 'linewidth',2)
% legend('RB')
%  figure(5)
%  plot(t*365, x3,'b', 'linewidth',2)
%  xlim([0,365])
%  xlabel('time (days')
%  legend('SA')
% figure(6)
% plot(t, x4,'b--', 'linewidth',2)
% legend('SB')
%  figure(7)
%  plot(t, x1,'r',t, x3,'g',t, x5,'b', 'linewidth',2)
%  legend('IA','SA','RA')
%  xlabel('time')
%  ylabel('individuals')
% figure(8)
% plot(t, x2,'r',t, x4,'g',t, x6,'b', 'linewidth',2)
% legend('IB','SB','RB')
% xlabel('time')
% ylabel('individuals')


%  xlim([0,50]);


%% Model equations

function [dX]=cholera(t,y)  %%phebl=patchy_ebola
% Parameters:
% Patch 1
mu1 = 1E-4; 
beta1I=2.64E-6;
beta1W=1.21E-4;
m1=.1; % S,R movement
delta1=5E-4;
gamma1=.25;
n1=0.001; % I movement
xi1=7.56E-3;
nu1=xi1; % decay of pathogen
rho1=3E-5; % pathogen transport rate
% Patch 2
mu2=mu1; 
beta2I=beta1I;
beta2W=beta1W;
m2=m1;
delta2=delta1;
gamma2=.25;
n2=n1;
xi2=xi1;
nu2=xi1;
rho2=3E-6;
%Variables: y1=S1, y2=I1, y3=R1, y4=W1, y5=S2, y6=I2, y7=R2, y8=W2
N1 = y(1) + y(2) + y(3);
N2 = y(5) + y(6) + y(7);

dX=zeros(8,1);
% Susceptible 1
dX(1) = mu1*N1 - beta1I*y(1)*y(2) - beta1W*y(1)*y(4) - mu1*y(1) - m1*y(1) + m2*y(5);
% Infectious 1
dX(2) = beta1I*y(1)*y(2) + beta1W*y(1)*y(4) - (gamma1 + mu1 + delta1)*y(2) - n1*y(2) + n2*y(6);
% Recovered 1
dX(3) = gamma1*y(2) - mu1*y(3) - m1*y(3) + m2*y(7);
% Water 1
dX(4) = xi1*y(2) - nu1*y(4) - rho1*y(4) + rho2*y(8);
% Susceptible 2
dX(5) = mu2*N2 - beta2I*y(5)*y(6) - beta2W*y(5)*y(8) - mu2*y(5) + m1*y(1) - m2*y(5);
% Infectious 2
dX(6) = beta2I*y(5)*y(6) + beta2W*y(5)*y(8) - (gamma2 + mu2 + delta2)*y(6) + n1*y(2) - n2*y(6);
% Recovered 2
dX(7) = gamma2*y(6) - mu2*y(7) + m1*y(3) - m2*y(7);
% Water 2
dX(8) = xi2*y(6) - nu2*y(8) + rho1*y(4) - rho2*y(8);
end