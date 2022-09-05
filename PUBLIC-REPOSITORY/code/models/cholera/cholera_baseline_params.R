# define parameters --------------------------------------------------------
params_cholera <- data.frame(
  ### ODE parameters
  mu1 = 0, mu2 = 0, # natural birth/death rate 1E-4
  beta_I1 = 2.64E-6, beta_I2 = 2.64E-6, # transmission rate from people (consider settings to 0)
  beta_W1 = 1.01E-5, beta_W2 = 1.01E-5, # transmission rate from water (consider increasing an order of magnitude)
  v1 = 0, v2 = 0, # vaccination rate
  u1 = 0, u2 = 0, # reduction in transmission due to sanitation
  m1 = 5e-4, m2 = 5e-4, # movement rate (non-infected)
  n1 = 0, n2 = 0, # movement rate (infected)
  gamma1 = 0.25, gamma2 = 0.25, # recovery rate
  delta1 = 5E-4, delta2 = 5E-4, # disease induced mortality
  xi1 = 7.56E-3, xi2 = 7.56E-3, # pathogen survival rate in water
  nu1 = 7.56E-3, nu2 = 7.56E-3, # pathogen clearance rate in water
  rho1 = 0.0125, rho2 = 0.0125, # pathogen movement rate in water (consider decreasing by half)
  ### optimal control parameters
  b1 = 1, b2 = 1, # cost of cases
  C1 = 0.125, C2 = 0.125, # cost of vaccinations
  epsilon1 = 10000, epsilon2 = 10000, # non-linearity for vacc
  D1 = 0.0125, D2 = 0.0125, # cost of sanitation
  eta1 = 100, eta2 = 100, # non-linearity for sanitation
  tol = 0.001, # optimization tolerance
  control_type = "unique",
  ### optimal control bounds
  v1_min = 0, v1_max = 0.015,
  v2_min = 0, v2_max = 0.015,
  u1_min = 0, u1_max = 0.4,
  u2_min = 0, u2_max = 0.4
)

# define time series (units of days)
times_cholera <- seq(0, 200, 0.05)

# initial control guesses
guess_v1 <- rep(0, length(times_cholera))
guess_v2 <- rep(0, length(times_cholera))
guess_u1 <- rep(0, length(times_cholera))
guess_u2 <- rep(0, length(times_cholera))


# define initial conditions (ICs)-----------------------------------------------
# run uncontrolled outbreak (beginning with ) for response time days,
# use the states on this day

response_time_cholera <- 60 # define time of outbreak response in days

# use initially uncontrolled outbreak to determine initial conditions
IC_cholera_uncontrol <- c(
  S1 = 100000 - 100, S2 = 100000, # consider doubling population of Patch 1
  I1 = 100, I2 = 0,
  R1 = 0, R2 = 0,
  W1 = 0, W2 = 0
)