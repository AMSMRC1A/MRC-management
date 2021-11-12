# define ODE parameters
params <- data.frame(
  mu1 = 0, mu2 = 0, # natural birth/death rate 1E-4
  beta_I1 = 2.64E-5, beta_I2 = 2.64E-5, # transmission rate from people (consider setting to 0)
  beta_W1 = 1.01E-4, beta_W2 = 1.01E-4, # transmission rate from water (consider increasing by an order of magnitude)
  v1 = 0, v2 = 0, # vaccination rate
  m1 = 0, m2 = 0, # movement rate (non-infected)
  n1 = 0, n2 = 0, # movement rate (infected)
  gamma1 = 0.25, gamma2 = 0.25, # recovery rate
  delta1 = 5E-4, delta2 = 5E-4, # disease induced mortality
  xi1 = 7.56E-3, xi2 = 7.56E-3, # pathogen survival rate in water
  nu1 = 7.56E-3, nu2 = 7.56E-3, # pathogen clearance rate in water
  rho1 = 0.025, rho2 = 0.025 # pathogen movement rate in water (consider decreasing by half)
)

# bounds of optimal control (OC) parameters
bounds <- list(MaxV1 = 0.015, MaxV2 = 0.015)
max_params <- params
max_params$v1 <- bounds$M1
max_params$v2 <- bounds$M2

# define time series (units of days)
times <- seq(0, 200, 0.05)

# define initial conditions (ICs)

# # set initial conditions explicitly
# IC <- c(S1 = 10000-100, S2 = 10000-10, # consider doubling population of Patch 1
#         I1 = 100,    I2 = 10,
#         R1 = 0,      R2 = 0,
#         W1 = 100,    W2 = 100)#,
# # c1 = 0,      c2 = 0, # added c to count new cases
#  V1 = 0,      V2 = 0) # added V to count vaccinations

# use initially uncontrolled outbreak to determine initial conditions
response_time <- 50 # define time of outbreak response in days
# See "Cholera_runODE.R" for details on how this time series was calculated
uncontrolled_epidemic <- read.csv("analysis/initialOutbreak_noControl_forIC.csv")
IC <- as.double(uncontrolled_epidemic[uncontrolled_epidemic$time == response_time, -1])
names(IC) <- names(uncontrolled_epidemic[, -1])
