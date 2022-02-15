N1<-1e5
N2<-1e5

##Set-up parameters
parm<-list(mu1=5.5e-5, #Change average lifespan to 50 years from 27
        alpha1=.1,
        gammaI1=1/15,
        gammaH1=.028,
        phi1=.236,
        deltaI1=.024,
        deltaH1=.01,
        xi1=0.222,
        #m1=0,
        m1=0, #movement
        n1=0,
        #Patch 2
        mu2=5.5e-5, #Change average lifespan to 50 years from 27
        alpha2=.1,
        gammaI2=1/15,
        gammaH2=.028,
        phi2=.236,
        deltaI2=.024,
        deltaH2=.01,
        xi2=0.222,
        #m2=0,
        m2=0.05,  #movement
        n2=0,
        b1=1,
        b2=1,
        C1=.01,
        C2=.01,
        epsilon1=5e7,  #5e7
        epsilon2=5e7   #5e7 for 2  #1e4 for 1
)

#Calculate betas based on R0
R0 <- 1.7
p <- 10  #scale betaI for transmission from D to get betaD

#patch 1
betaI1 <- .006/N1*1.8 #Burton
#betaI1 <-  with(parm,R0/(N1*(alpha1/(alpha1+mu1))*(1/(gammaI1+phi1+deltaI1+mu1))+p*N1*(alpha1/(alpha1+mu1))*(deltaI1/(gammaI1+phi1+deltaI1+mu1))*(1/(xi1)))/1)#
betaD1 <- 3.3/N1*1.8 #Burton
#betaD1 <- p*betaI1


#patch 2
betaI2 <- .006/N1*1.8 #Burton
#betaI2 <-  with(parm,R0/(N2*(alpha2/(alpha2+mu2))*(1/(gammaI2+phi2+deltaI2+mu2))+p*N2*(alpha2/(alpha2+mu2))*(deltaI2/(gammaI2+phi2+deltaI2+mu2))*(1/(xi2)))/1)
betaD2 <- 3.3/N1*1.8 #Burton
#betaD2 <- p*betaI2


parm2<-parm
parm<-c(parm,betaI1=betaI1,betaI2=betaI2,betaD1=betaD1,betaD2=betaD2,N1=N1,N2=N2)


#Initial conditions
Y=c(S1=N1-700,E1=400,I1=300,H1=0,D1=0,R1=0,S2=N2,E2=0,I2=0,H2=0,D2=0,R2=0)
rm(betaI1,betaI2,betaD1,betaD2,N1,N2,R0,p)
