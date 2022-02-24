# define parameters --------------------------------------------------------
# population sizes
N1<-1e5
N2<-1e5

params_ebola <- data.frame(
  mu1=5.5e-5, #Change average lifespan to 50 years from 27
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
betaD1 <- 3.3/N1*1.8 #Burton
#patch 2
betaI2 <- .006/N1*1.8 #Burton
betaD2 <- 3.3/N1*1.8 #Burton

params_ebola<-c(params_ebola,betaI1=betaI1,betaI2=betaI2,betaD1=betaD1,betaD2=betaD2,N1=N1,N2=N2)

# define time series (units of days)
times_ebola <- seq(0, 100, 0.05)

#initial control guesses
guess_v1 <- rep(0, length(times_ebola))
guess_v2 <- rep(0, length(times_ebola))
guess_u1 <- rep(0, length(times_ebola))
guess_u2 <- rep(0, length(times_ebola))


# define initial conditions (ICs)-----------------------------------------------
IC_ebola <- c(
  S1=N1-700, S2=N2,
  E1=400,    E2=0,
  I1=300,    I2=0,
  H1=0,      H2=0,
  D1=0,      D2=0,
  R1=0,      R2=0
)

rm(betaI1,betaI2,betaD1,betaD2,N1,N2,R0,p)
