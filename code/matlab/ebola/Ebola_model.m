function dy = Ebola_model(t, y, par)

% pull out parameters

betaI1 = par.betaI1;
betaI2 = par.betaI2;
betaD1 = par.betaD1;
betaD2 = par.betaD2;
alpha1 = par.alpha1;
alpha2 = par.alpha2;
gammaI1 = par.gammaI1;
gammaI2 = par.gammaI2;
gammaH1 = par.gammaH1;
gammaH2 = par.gammaH2;
phi1 = par.phi1;
phi2 = par.phi2;
deltaI1 = par.deltaI1;
deltaI2 = par.deltaI2;
deltaH1 = par.deltaH1;
deltaH2 = par.deltaH2;
xi1 = par.xi1;
xi2 = par.xi2;
mu1 = par.mu1;
mu2 = par.mu2;
v1 = par.v1;
v2 = par.v2;
m1 = par.m1;
m2 = par.m2;
n1 = par.n1;
n2 = par.n2;

% define states

S1 = y(1);
E1 = y(2);
I1 = y(3);
H1 = y(4);
D1 = y(5);
R1 = y(6);

S2 = y(7);
E2 = y(8);
I2 = y(9);
H2 = y(10);
D2 = y(11);
R2 = y(12);

N1 = S1 + E1 + I1 + H1 + R1;
N2 = S2 + E2 + I2 + H2 + R2;


% define ODEs

% cumulative cases
dci = betaI1*S1*I1 + betaD1*S1*D1;

% patch 1
dS1 = mu1*N1 - betaI1*S1*I1 - betaD1*S1*D1 - mu1*S1 - v1*S1 - m1*S1 + m2*S2 ;
dE1 = betaI1*S1*I1 + betaD1*S1*D1 - mu1*E1 - alpha1*E1 - m1*E1 + m2*E2 ;
dI1 = alpha1*E1 - (mu1 + gammaI1 + phi1 + deltaI1)*I1 - n1*I1 + n2*I2 ;
dH1 = phi1*I1 - (gammaH1 + deltaH1 + mu1)*H1 ;
dD1 = deltaI1*I1 - xi1*D1 ;
dR1 = gammaI1*I1 + gammaH1*H1 + v1*S1 - mu1*R1 - m1*R1 + m2*R2 ;

% patch 2
dS2 = mu2*N2 - betaI2*S2*I2 - betaD2*S2*D2 - mu2*S2 - v2*S2 + m1*S1 - m2*S2 ;
dE2 = betaI2*S2*I2 + betaD2*S2*D2 - mu2*E2 - alpha2*E2 + m1*E1 - m2*E2 ;
dI2 = alpha2*E2 - (mu2 + gammaI2 + phi2 + deltaI2)*I2 + n1*I1 - n2*I2 ;
dH2 = phi2*I2 - (gammaH2 + deltaH2 + mu2)*H2 ;
dD2 = deltaI2*I2 - xi2*D2 ;
dR2 = gammaI2*I2 + gammaH2*H2 + v2*S2 - mu2*R2 + m1*R1 - m2*R2 ;

% combine cases
dy = [dS1, dE1, dI1, dH1, dD1, dR1, dS2, dE2, dI2, dH2, dD2, dR2, dci]'; 


end