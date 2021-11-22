

##Set-up parameters
p<-list(mu1=1e-4,
        betaI1=10e-9,
        betaD1=10e-6,
        alpha1=.1,
        gammaI1=.02,
        gammaH1=.028,
        phi1=.41,
        deltaI1=.024,
        deltaH1=.01,
        xi1=0.222,
        v1=0,
        m1=0,
        n1=0,
        #Patch 2
        mu2=1e-4,
        betaI2=10e-9,
        betaD2=10e-6,
        alpha2=.1,
        gammaI2=.02,
        gammaH2=.028,
        phi2=.41,
        deltaI2=.024,
        deltaH2=.01,
        xi2=0.222,
        v2=0,
        m2=0,
        n2=0
)

#Initial conditions
Y=c(S1=(1e6)-1,E1=0,I1=1,H1=0,D1=0,R1=0,S2=(1e6),E2=0,I2=0,H2=0,D2=0,R2=0)

