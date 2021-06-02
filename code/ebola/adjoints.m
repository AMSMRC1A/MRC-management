function [ dlambda ] = adjoints(t,lambda,tvec,x,v1,v2,par)

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
m1 = par.m1;
m2 = par.m2;
n1 = par.n1;
n2 = par.n2;
b1 = par.b1;
b2 = par.b2;
C1 = par.C1;
C2 = par.C2;

% define states
S1 = x(1);
E1 = x(2);
I1 = x(3);
H1 = x(4);
D1 = x(5);
R1 = x(6);

S2 = x(7);
E2 = x(8);
I2 = x(9);
H2 = x(10);
D2 = x(11);
R2 = x(12);

% define states
lambda1 = lambda(1);
lambda2 = lambda(2);
lambda3 = lambda(3);
lambda4 = lambda(4);
lambda5 = lambda(5);
lambda6 = lambda(6);

lambda7 = lambda(7);
lambda8 = lambda(8);
lambda9 = lambda(9);
lambda10 = lambda(10);
lambda11 = lambda(11);
lambda12 = lambda(12);

% define ODEs
x=interp1(tvec,x,t);
v1=pchip(tvec,v1,t);
v2=pchip(tvec,v2,t);
dlambda=zeros(12,1);

% patch 1
dlambda(1) = -( b1*(betaI1*I1 + betaD1*D1) + C1*v1 + lambda1*(mu1 -betaI1*I1 - betaD1*D1 - mu1 - v1 - m1)  + lambda2*(betaI1*I1 + betaD1*D1) + lambda6*v1 + lambda7*m1 );
dlambda(2) =-( C1*v1  + lambda1*mu1 + lambda2*(-mu1 - alpha1 - m1) + lambda3*alpha1 + lambda8*m1 );
dlambda(3) = -( b1*(betaI1*S1) + lambda1*(mu1 - betaI1*S1) + lambda2*(betaI1*S1) + lambda3*(-(mu1 + gammaI1 + phi1 + deltaI1 + n1))  + lambda4*phi1 + lambda5*deltaI1 + lambda6*gammaI1 + lambda9*(n1) );
dlambda(4) = -(lambda1*(mu1) + lambda4*(-(mu1 + gammaH1 + deltaH1)) + lambda6*(gammaH1));
dlambda(5) =-( b1*(betaD1*S1) - lambda1*betaD1*S1 + lambda2*betaD1*S1 - lambda5*xi1);
dlambda(6) = -( lambda1*mu1 + lambda6*(-mu1 - m1) + lambda12*m1);

% patch 2
dlambda(7) = -( b2*(betaI2*I2 + betaD2*D2) + C2*v2 + lambda7*(mu2 -betaI2*I2 - betaD2*D2 - mu2 - v2 - m2)  + lambda8*(betaI2*I2 + betaD2*D2) + lambda12*v2 + lambda1*m2);
dlambda(8) =-( C2*v2  + lambda7*mu2 + lambda8*(-mu2 - alpha2 - m2) + lambda9*alpha2 + lambda2*m2);
dlambda(9) =- (b2*betaI2*S2 + lambda7*(mu2 - betaI2*S2) + lambda8*(betaI2*S2) + lambda9*(-(mu2 + gammaI2 + phi2 + deltaI2 + n2)) + lambda10*phi2 + lambda11*deltaI2 + lambda12*gammaI2 + lambda3*n2); 
dlambda(10) =-(lambda7*mu2 + lambda10*(-(mu2 + gammaH2 + deltaH2)) + lambda12*gammaH2);
dlambda(11) = -( b2*betaD2*S2 - lambda7*betaD2*S2 + lambda8*betaD2*S2 - lambda11*xi2 );
dlambda(12)= -( lambda7*mu2 + lambda12*(-mu2 - m2) + lambda6*m2);


end

