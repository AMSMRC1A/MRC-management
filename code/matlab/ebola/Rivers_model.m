function dy = Rivers_model(t, y, p)

% define states

S1 = y(1);
E1 = y(2);
I1 = y(3);
H1 = y(4);
D1 = y(5);
R1 = y(6);

%S2 = y(7);
%E2 = y(8);
%I2 = y(9);
%H2 = y(10);
%D2 = y(11);
%R2 = y(12);

N1 = S1 + E1 + I1 + H1 + D1 + R1;
%N2 = S2 + E2 + I2 + H2 + D1 + R2;


% define ODEs

% cumulative cases
dci = p.betaI*S1*I1/N1 + p.betaF*S1*D1/N1 + p.betaH*S1*H1/N1;

% patch 1
dS1 = -dci;
dE1 = dci - p.alpha*E1;
dI1 = p.alpha*E1 - (p.gammaH*p.theta1 + p.gammaI*(1-p.theta1)*(1-p.delta1) + p.gammaD*(1-p.theta1)*p.delta1)*I1;
dH1 = p.gammaH*p.theta1*I1 - (p.gammaDH*p.delta2 + p.gammaIH*(1-p.delta2))*H1;
dD1 = p.gammaD*(1-p.theta1)*p.delta1*I1 + p.gammaDH*p.delta2*H1 - p.gammaD*D1;
dR1 = p.gammaI*(1-p.theta1)*(1-p.delta1)*I1 + p.gammaIH*(1-p.delta2)*H1 + p.gammaD*D1;

% patch 2
%dS2 = mu2*N2 - betaI2*S2*I2 - betaD2*S2*D2 - mu2*S2 - v2*S2 + m1*S1 - m2*S2 ;
%dE2 = betaI2*S2*I2 + betaD2*S2*D2 - mu2*E2 - alpha2*E2 + m1*E1 - m2*E2 ;
%dI2 = alpha2*E2 - (mu2 + gammaI2 + phi2 + deltaI2)*I2 + n1*I1 - n2*I2 ;
%dH2 = phi2*I2 - (gammaH2 + deltaH2 + mu2)*H2 ;
%dD2 = deltaI2*I2 - xi2*D2 ;
%dR2 = gammaI2*I2 + gammaH2*H2 + v2*S2 - mu2*R2 + m1*R1 - m2*R2 ;

% combine cases
%dy = [dS1, dE1, dI1, dH1, dD1, dR1, dS2, dE2, dI2, dH2, dD2, dR2, dci]'; 

dy = [dS1, dE1, dI1, dH1, dD1, dR1, dci]'; 

end