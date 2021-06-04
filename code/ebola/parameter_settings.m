

%% Initial Conditions

% patch 1
N1 = 6e5;
S1 = N1 - 100;
E1 = 50;
I1 = 50;
H1 = 0;
D1 = 0;
R1 = 0;
ci = 0;

% patch 2
N2 = 6e5;
S2 = N2;
E2 = 0;
I2 = 0;
H2 = 0;
D2 = 0;
R2 = 0;

% combine into single vector
y0 = [S1, E1, I1, H1, D1, R1, S2, E2, I2, H2, D2, R2];



%% All parameters for two-patch Ebola model 

alpha1 = 0.1;
alpha2 = 0.1;
gammaI1 = 1/15; %originally 0.02 (Burton)
gammaI2 = 1/15; %originally 0.02 (Burton)
gammaH1 = 0.028;
gammaH2 = 0.028;
phi1 = 0.236; %average from Burton
phi2 = 0.236;
deltaI1 = 0.024;
deltaI2 = 0.024;
deltaH1 = 0.01;
deltaH2 = 0.01;
xi1 = 0.222;
xi2 = 0.222;
mu1 = 1e-4;
mu2 = 1e-4;
%v1 = 0;
%v2 = 0;
m1 = 0.005;
m2 = 0.005;
n1 = 0;
n2 = 0;
b1 = 1;
b2 = 1;
C1 = 0.5;
C2 = 0.5;
epsilon1 = 0.5;
epsilon2 = 0.5;

% calculate betas based on R0
R0 = 1.7;
p = 10; % scale betaI for transmission from D to get betaD % from Julie and Lauren's paper

% patch 1
%betaI1 = 1e-9; % Burton
betaI1 = R0/(N1*(alpha1/(alpha1+mu1))*(1/(gammaI1+phi1+deltaI1+mu1))+p*N1*(alpha1/(alpha1+mu1))*(deltaI1/(gammaI1+phi1+deltaI1+mu1))*(1/(xi1)))/1;
%betaD1 = 5.5e-7; % Burton
betaD1 = p*betaI1;

% patch 2
%betaI2 = 1e-9;
betaI2 = betaI1;
%betaD2 = 5e-7;
betaD2 = betaD1;


% Add parameters to structure 'par'
par.betaI1 = betaI1;
par.betaI2 = betaI2;
par.betaD1 = betaD1;
par.betaD2 = betaD2;
par.alpha1 = alpha1;
par.alpha2 = alpha2;
par.gammaI1 = gammaI1;
par.gammaI2 = gammaI2;
par.gammaH1 = gammaH1;
par.gammaH2 = gammaH2;
par.phi1 = phi1;
par.phi2 = phi2;
par.deltaI1 = deltaI1;
par.deltaI2 = deltaI2;
par.deltaH1 = deltaH1;
par.deltaH2 = deltaH2;
par.xi1 = xi1;
par.xi2 = xi2;
par.mu1 = mu1;
par.mu2 = mu2;
%par.v1 = v1;
%par.v2 = v2;
par.m1 = m1;
par.m2 = m2;
par.n1 = n1;
par.n2 = n2;
par.b1 = b1;
par.b2 = b2;
par.C1 = C1;
par.C2 = C2;
par.epsilon1 = epsilon1;
par.epsilon2 = epsilon2;









