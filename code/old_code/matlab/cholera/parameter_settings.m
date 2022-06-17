%% Parameter values 
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

%% Initial conditions
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
