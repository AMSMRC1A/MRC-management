# Cholera optimal control "difference optimization" function

# Load necessary packages-------------------------------------------------------
# library(nloptr)
library(foreach)
library(pracma)
library(optimization)
library(tidyverse)

# load files--------------------------------------------------------------------
source("CholeraSIRW_ODE.R")
source("Cholera_params.R")
source("Cholera_OCfunc.R")
source("Cholera_adjoints.R")
source("Cholera_analysisfunc.R")

# Optimization set up-----------------------------------------------------------
# Inputs:
# - parameter values (initial guesses)
# - parameters that can be varied
# - parameter ranges/bounds

# Objective function:
# - relative difference in cost between uniform and unique controls
#   get.rel.cost function

# Output:
# - parameter set which maximizes relative cost difference
# - value of the relative cost difference



# Define baseline parameters----------------------------------------------------
# Baseline values derived from "Cholera_params.R"
# These are: params, IC, and times

# Define variable parameters and their ranges-----------------------------------
# Initial conditions
# ** Can't change these with the way apply_oc is set up **
# Total population sizes


# Initial no. infected individuals

# Initial amount of pathogen in water

# Transmission parameters
# betaI's

# betaW's

# Optimal control settings------------------------------------------------------
# initial guesses - controls
guess_v1 <- rep(0, length(times))
guess_v2 <- rep(0, length(times))

# Cost of cases
tol = 0.01 # tolerance parameter for optimization
oc_params <- c(b1 = 1, b2 = 1, # cost of cases
               C1 = 0.125, C2 = 0.125,  # cost of vaccinations
               epsilon1  = 10000, epsilon2 = 10000) # non-linearity

# Set up parameter tables-------------------------------------------------------
## ** The following might be circular **

temp_params <- left_join(as_tibble(as.list(params)),
                         as_tibble(as.list(oc_params)),
                         by = character()) %>% 
  # select which parameters you want to vary
  select(m1, m2, beta_I1, beta_I2, beta_W1, beta_W2, rho1, rho2, b2, C2)
temp_params$m1 <- 0.025
temp_params$m2 <- 0.025

variable.parameters <- tibble(
  var_name = names(temp_params),
  # initial values are determined from our original parameterization
  value = as.double(temp_params[1,]),
  # bounds are stuff I made up
  lower_bounds = c(0, 0, temp_params$beta_I1*1E-1, temp_params$beta_I2*1E-1,
                   temp_params$beta_W1*1E-1, temp_params$beta_W2*1E-1, 0, 0,
                   0.1, temp_params$C2*1E-1),
  upper_bounds = c(1, 1, temp_params$beta_I1*1E1, temp_params$beta_I2*1E1,
                   temp_params$beta_W1*1E1, temp_params$beta_W2*1E1, 0.05, 0.05,
                   10,temp_params$C2*1E1)
)

# Turn the variable.parameters data frame into something that can be put into the "apply_oc" function
num_vars <- dim(variable.parameters)[1]
test_params <- tibble(!!variable.parameters$var_name[1] := !!variable.parameters$value[1])
for (i in 2:num_vars) {
  test_params <- add_column(test_params,
                            !!variable.parameters$var_name[i] := !!variable.parameters$value[i])
}
test_params <- as.data.frame(test_params)


# Functions for optimization: --------------------------------------------------
# get.rel.costs: get the relative cost difference between using uniform and unique controls
get.rel.cost <- function(test_params, values) {
  # Get test params
  # test_params <- helper_df_func(variable.parameters)
  # test_params[1,] <- values
  test_params <- expand_grid(test_params, control_type = c("unique", "uniform"))
  
  # Run optimal control calculations across test_params dataframe
  vary_params <- foreach(
    i = 1:nrow(test_params),
    .packages = c("deSolve", "tidyverse", "pracma"),
    # explicitly give 'foreach' the functions and data it needs
    .export = c(
      "apply_oc", "guess_v1", "guess_v2", "IC",
      "bounds", "chol", "adj", "times",
      "params", "oc_params", "run_oc", "tol",
      "oc_optim", "calc_opt_v", "norm_oc", "calc_j",
      "eval_j_integrand"
    ),
    .combine = rbind
  ) %dopar% {
    apply_oc(
      change_params = test_params[i, ],
      guess_v1 = guess_v1, guess_v2 = guess_v2,
      init_x = IC, bounds = bounds,
      ode_fn = chol, adj_fn = adj,
      times = times, params = c(params, oc_params),
      tol = tol,
      control_type = test_params[i, "control_type"],
      return_type = "j"
    )
  }
  # reformat to a data-frame
  reformat_df <- reformat_mult_params_output(vary_params, test_params)
  # compute relative cost
  rel.cost <- reformat_df$j %>%
    mutate(j_change = 100 * (j - j[control_type == "uniform"]) / j[control_type == "uniform"]) %>%
    filter(control_type == "unique") %>%
    select(j_change)
  return(rel.cost$j_change)
}

# Optimization function
# Same as above, but sets test_params = test_params
optim_func <- function(x) {get.rel.cost(test_params,x)}

# Optimization problem----------------------------------------------------------
# Using 'optim' package
res_optim <- optim(
  par = variable.parameters$value,
  fn = optim_func,
  method = "L-BFGS-B",
  lower = variable.parameters$lower_bounds,
  upper = variable.parameters$upper_bounds,
  control=list(trace=6)
)
res_optim_SANN <- optim( # simulated annealing. doesn't work with constraints
  par = variable.parameters$value,
  fn = optim_func,
  method = "SANN",
  #lower = variable.parameters$lower_bounds,
  #upper = variable.parameters$upper_bounds,
  control=list(trace=6)
)
# Didn't really work.

# Using 'nloptr' package
# opts <- list("algorithm"="NLOPT_GN_DIRECT_NOSCAL",
#              "xtol_rel"=1.0e-8)
# res_nloptr <- nloptr(
#   x0 = variable.parameters$value,
#   eval_f = optim_func,
#   lb = variable.parameters$lower_bounds,
#   ub = variable.parameters$upper_bounds,
#   opts = opts
# )
#

# Using 'optimization' package
res_optim_sa <- optim_sa(fun = optim_func,
                         start = variable.parameters$value,
                         trace = TRUE,
                         lower = variable.parameters$lower_bounds,
                         upper = variable.parameters$upper_bounds,
                         control = list(t0 = 500, nlimit = 50, r = 0.85,
                                        rf = c(1,1,1,1,1,1,1,1,1,1), ac_acc = 0.1, dyn_rf = TRUE))

## NOTES
# Only varying m1 and m2, the largest rel. cost diff. is 1%, which occurs when m1=1 and m2=0
