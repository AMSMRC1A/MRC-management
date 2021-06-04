N1<-6e4
N2<-6e4

##Set-up parameters
parm<-list(mu1=1e-4,
        alpha1=.1,
        gammaI1=.02,
        gammaH1=.028,
        phi1=.41,
        deltaI1=.024,
        deltaH1=.01,
        xi1=0.222,
        m1=0.005,
        n1=0,
        #Patch 2
        mu2=1e-4,
        alpha2=.1,
        gammaI2=.02,
        gammaH2=.028,
        phi2=.41,
        deltaI2=.024,
        deltaH2=.01,
        xi2=0.222,
        m2=0.005,
        n2=0,
        b1=1,
        b2=1,
        C1=.95,
        C2=.95,
        epsilon1=.01,
        epsilon2=.01
)

#Calculate betas based on R0
R0 <- 1.7
p <- 10  #scale betaI for transmission from D to get betaD

#patch 1
#betaI1 <- 1e-9 #Burton
betaI1 <-  with(parm,R0/(N1*(alpha1/(alpha1+mu1))*(1/(gammaI1+phi1+deltaI1+mu1))+p*N1*(alpha1/(alpha1+mu1))*(deltaI1/(gammaI1+phi1+deltaI1+mu1))*(1/(xi1)))/1)#
#betaD2 <- 5e-7 #Burton
betaD1 <- p*betaI1


#patch 2
#betaI2 <- 1e-9 #Burton
betaI2 <-  with(parm,R0/(N2*(alpha2/(alpha2+mu2))*(1/(gammaI2+phi2+deltaI2+mu2))+p*N2*(alpha2/(alpha2+mu2))*(deltaI2/(gammaI2+phi2+deltaI2+mu2))*(1/(xi2)))/1)
#betaD2 <- 5e-7 #Burton
betaD2 <- p*betaI2

parm<-c(parm,betaI1=betaI1,betaI2=betaI2,betaD1=betaD1,betaD2=betaD2)

#Initial conditions
Y=c(S1=N1-1000,E1=600,I1=400,H1=0,D1=0,R1=0,S2=N2,E2=0,I2=0,H2=0,D2=0,R2=0)
rm(betaI1,betaI2,betaD1,betaD2,N1,N2,R0,p)
