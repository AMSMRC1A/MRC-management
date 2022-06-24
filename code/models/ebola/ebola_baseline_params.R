# define parameters --------------------------------------------------------
# population sizes
N1<-1e5
N2<-1e5

params_ebola <- data.frame(
  mu1=5.5e-5, #Change average lifespan to 50 years from 27
  alpha1=.1,
  gammaI1=1/10, #1/15,
  gammaH1=.154,
  phi1=.236,
  deltaI1=.024,
  deltaH1=.01,
  xi1=0.222,
  #m1=0,
  m1=5e-4, #movement
  n1=0,
  #Patch 2
  mu2=5.5e-5, #Change average lifespan to 50 years from 27
  alpha2=.1,
  gammaI2=1/10, #1/15,
  gammaH2=.154,
  phi2=.236,
  deltaI2=.024,
  deltaH2=.01,
  xi2=0.222,
  #m2=0,
  m2=5e-4,  #movement
  n2=0,
  b1=1,
  b2=1,
  Cv1=.01,
  Cv2=.01,
  epsilonV1=5e4,  #5e7
  epsilonV2=5e4, #5e7 for 2  #1e4 for 1  
  Cu1=0.1, #.001,
  Cu2=0.1, #.001,
  epsilonU1=5E1, #5e4, 5e5?
  epsilonU2=5E1, #5e4, 5e5?
  # for baseline control_type
  control_type = "uniform",
  tol= 0.01, # optimization tolerance
  v1 = 0, v2 = 0, # baseline vaccination rate
  u1 = 0, u2 = 0, # baseline sanitation rate
  ### optimal control bounds
  v1_min = 0, v1_max = 0.015, 
  v2_min = 0, v2_max = 0.015,
  u1_min = 0, u1_max = 0.5, #phi1 
  u2_min = 0, u2_max = 0.5  #phi2
)

#Calculate betas based on R0
R0 <- 1.7
p <- 10  #scale betaI for transmission from D to get betaD

#patch1
betaI1 <- with(params_ebola, R0/(N1*(alpha1/(alpha1+mu1))*(1/(gammaI1+phi1+deltaI1+mu1)) + p*N1*(alpha1/(alpha1+mu1))*(deltaI1/(gammaI1+phi1+deltaI1+mu1))*(1/xi1)))
betaD1 <- p*betaI1
#patch2
betaI2 <- betaI1
betaD2 <- betaD1

params_ebola<-c(params_ebola,betaI1=betaI1,betaI2=betaI2,betaD1=betaD1,betaD2=betaD2,N1=N1,N2=N2)

# define time series (units of days)
times_ebola <- seq(0, 200, 0.05)

#initial control guesses
guess_v1 <- rep(0, length(times_ebola))
guess_v2 <- rep(0, length(times_ebola))
guess_u1 <- rep(0, length(times_ebola))
guess_u2 <- rep(0, length(times_ebola))


# define initial conditions (ICs)-----------------------------------------------
# IC_ebola <- c(
#   S1=N1-700, 
#   E1=400,    
#   I1=300,    
#   H1=0,      
#   D1=0,      
#   R1=0,      
#   S2=N2,
#   E2=0,
#   I2=0,
#   H2=0,
#   D2=0,
#   R2=0
# )
# 
# rm(betaI1,betaI2,betaD1,betaD2,N1,N2,R0,p)


# define initial conditions (ICs)-----------------------------------------------
# run uncontrolled outbreak (beginning with ) for response time days, 
# use the states on this day

response_time <- 150 # define time of outbreak response in days
# use initially uncontrolled outbreak to determine initial conditions
IC_init <- c(
  S1 = N1 - 10,
  E1 = 0,
  I1 = 10,
  H1 = 0,
  D1 = 0,
  R1 = 0,
  S2 = N2,
  E2 = 0,
  I2 = 0,
  H2 = 0,
  D2 = 0,
  R2 = 0
)

source("models/ebola/ebola_functions.R")
# solve ODE
uncontrolled <- ode(y = IC_init,
                    times = times_ebola,
                    func = ode_ebola,
                    parms = params_ebola)
# set IC based on response_time
IC_ebola <- as.double(uncontrolled[uncontrolled[,"time"] == response_time, -1])
names(IC_ebola) <- colnames(uncontrolled[, -1])

# as.data.frame(uncontrolled) %>%
#   reshape2::melt(id = c("time")) %>%
#   mutate(patch = substr(variable, 2,2), 
#          comp = substr(variable,1,1)) %>%
#   ggplot(aes(x = time, y =value,color = patch)) +
#   geom_line() +
#   facet_wrap(vars(comp), scales = "free")

