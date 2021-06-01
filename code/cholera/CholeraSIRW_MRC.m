%%%%%%%%%%% ODE solution for Cholera SIR two-patch model;%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  May 31st, 2021 %%%%%%%%%%%


% Parameter values 
mu1 = 1E-4; 
mu2=mu1; 
beta1I=2.64E-5;
beta2I=beta1I;
beta1W= 1.21E-4;
beta2W=beta1W;
m1=.1;
m2=m1;
delta1=5E-4;
delta2=delta1;
n1=0;
n2=n1;
gamma1=.25;
gamma2=.25;
xi1=7.56E-3;
xi2=xi1;
nu1=xi1;
nu2=xi1;
rho1=0;
rho2=.1;


%Reproduction number (ODE with constant transmission rates)
% F=[betAbar*SstA0/NA,0; 0,betBbar*SstB0/NB];
% V=[muA+gamaA+phiIAB, -phiIBA;-phiIAB, muB+gamaB+phiIBA];
% M=F*V^(-1);
% e=eig(M)
% % R0A=;%%expression is too long
% % R0B=;  
% R0=max(e) 




% Timespan in days
tspan=0:.1:60;

% Initial conditions
x10=5000;%S1
x20=0;%I1
x30=0;%R1
x40=0;%W1
x50=5000;%S2
x60=10;%I2
x70=0;%R2
x80=2;%W2
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
plot( t, y, 'linewidth',2)
legend()
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
% mu1, mu2, beta1I, beta2I, beta1W, beta2W, m1, m2, delta1, delta2, n1, n2,
% gamma1, gamma2, xi1, xi2, nu1, nu2, rho1, rho2
mu1 = 1E-4; 
mu2=mu1; 
beta1I=2.64E-5;
beta2I=beta1I;
beta1W= 1.21E-4;
beta2W=beta1W;
m1=.1;
m2=m1;
delta1=5E-4;
delta2=delta1;
n1=0;
n2=n1;
gamma1=.25;
gamma2=.25;
xi1=7.56E-3;
xi2=xi1;
nu1=xi1;
nu2=xi1;
rho1=0;
rho2=.1;
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