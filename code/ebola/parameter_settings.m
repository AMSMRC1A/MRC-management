

%% Initial Conditions

% patch 1
N1 = 6e6;
S1 = N1 - 100;
E1 = 0;
I1 = 100;
H1 = 0;
D1 = 0;
R1 = 0;
ci = 0;

% patch 2
N2 = 0;
S2 = 0;
E2 = 0;
I2 = 0;
H2 = 0;
D2 = 0;
R2 = 0;

% combine into single vector
y0 = [S1, E1, I1, H1, D1, R1, S2, E2, I2, H2, D2, R2, ci];



%% All parameters for two-patch Ebola model 

alpha1 = 0.1;
alpha2 = 0.1;
gammaI1 = 0.1; %originally 0.02 (Burton)
gammaI2 = 0.1; %originally 0.02 (Burton)
gammaH1 = 0.028;
gammaH2 = 0.028;
phi1 = 0.236;
phi2 = 0.236;
deltaI1 = 0.024;
deltaI2 = 0.024;
deltaH1 = 0.01;
deltaH2 = 0.01;
xi1 = 0.222;
xi2 = 0.222;
mu1 = 1e-4;
mu2 = 1e-4;
v1 = 0;
v2 = 0;
m1 = 0;
m2 = 0;
n1 = 0;
n2 = 0;

% calculate betas based on R0
R0 = 1.7;
p = 10; % scale betaI for transmission from D to get betaD

% patch 1
%betaI1 = 1e-9; % Burton
betaI1 = R0/(N1*(alpha1/(alpha1+mu1))*(1/(gammaI1+phi1+deltaI1+mu1))+p*N1*(alpha1/(alpha1+mu1))*(deltaI1/(gammaI1+phi1+deltaI1+mu1))*(1/(xi1)))/1;
%betaD1 = 5e-7; % Burton
betaD1 = p*betaI1;

% patch 2
betaI2 = 1e-9;
betaD2 = 5e-7;


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
par.v1 = v1;
par.v2 = v2;
par.m1 = m1;
par.m2 = m2;
par.n1 = n1;
par.n2 = n2;









