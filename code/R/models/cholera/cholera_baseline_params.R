# define parameters --------------------------------------------------------
params_cholera <- data.frame(
  ### ODE parameters
  mu1 = 0, mu2 = 0, # natural birth/death rate 1E-4
  beta_I1 = 2.64E-5, beta_I2 = 2.64E-5, # transmission rate from people (consider setting to 0)
  beta_W1 = 1.01E-4, beta_W2 = 1.01E-4, # transmission rate from water (consider increasing by an order of magnitude)
  v1 = 0, v2 = 0, # vaccination rate
  u1 = 0, u2 = 0, # reduction in transmission due to santation 
  m1 = 0, m2 = 0, # movement rate (non-infected)
  n1 = 0, n2 = 0, # movement rate (infected)
  gamma1 = 0.25, gamma2 = 0.25, # recovery rate
  delta1 = 5E-4, delta2 = 5E-4, # disease induced mortality
  xi1 = 7.56E-3, xi2 = 7.56E-3, # pathogen survival rate in water
  nu1 = 7.56E-3, nu2 = 7.56E-3, # pathogen clearance rate in water
  rho1 = 0.025, rho2 = 0.025, # pathogen movement rate in water (consider decreasing by half)
  ### optimal control parameters
  b1 = 1, b2 = 1, # cost of cases
  C1 = 0.125, C2 = 0.125, # cost of vaccinations
  epsilon1 = 10000, epsilon2 = 10000, # non-linearity for vacc
  D1 = 0.0125, D2 = 0.0125, # cost of sanitation
  eta1 = 1000, eta2 = 1000 # non-linearity for sanitation
)

# bounds of optimal control (OC) parameters
bounds_cholera <- list(
  V1_min = 0, V1_max = 0.015, 
  V2_min = 0, V2_max = 0.015, 
  U1_min = 0, U1_max = 0.4, 
  U2_min = 0, U2_max = 0.4)

# define time series (units of days)
times_cholera <- seq(0, 200, 0.05)

# define initial conditions (ICs)-----------------------------------------------
# run uncontrolled outbreak (beginning with ) for response time days, 
# use the states on this day

response_time <- 100 # define time of outbreak response in days
# use initially uncontrolled outbreak to determine initial conditions
IC_init <- c(
  S1 = 100000 - 1, S2 = 100000, # consider doubling population of Patch 1
  I1 = 1, I2 = 0,
  R1 = 0, R2 = 0,
  W1 = 0, W2 = 0
)

# solve ODE
uncontrolled <- ode(y = IC_init, 
           times = times_cholera, 
           func = ode_cholera, 
           parms = params_cholera, 
           v1_interp = v1_interp, 
           v2_interp = v2_interp, 
           u1_interp = u1_interp, 
           u2_interp = u2_interp)
# set IC based on response_time
IC_cholera <- as.double(uncontrolled[uncontrolled$time == response_time, -1])
names(IC_cholera) <- names(uncontrolled[, -1])