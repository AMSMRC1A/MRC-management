# Cholera optimal control "difference optimization" function

# Load necessary packages-------------------------------------------------------
library(nloptr)

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

# Helper function: get.costs----------------------------------------------------
# Take in a set of parameters and get the costs of using uniform or unique
# controls
get.rel.cost <- function(x) {
  # Expand input to include both control types
  test_params <- data.frame(
    m1 = c(x[1], x[1]),
    m2 = c(x[2], x[2]),
    control_type = c("unique", "uniform")
  ) # temporary. need to find a way to pass the variable names to this function

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

# Define baseline parameters----------------------------------------------------
# Baseline values derived from "Cholera_params.R"
# These are: params, IC, and times

# Define variable parameters and their ranges-----------------------------------
# Initial conditions
# Total population sizes


# Initial no. infected individuals

# Initial amount of pathogen in water

# Transmission parameters
# betaI's

# betaW's

# Optimal control parameters
# Cost of cases

# Cost of vaccination (relative to cases)

# Non-linear term (epsilon)

# Helper function to put variables togeth
variable.parameters <- data.frame(
  var_name = c("m1", "m2"),
  value = c(0.025, 0.025),
  lower_bounds = c(0, 0),
  upper_bounds = c(0.05, 0.05)
)

# Optimization problem----------------------------------------------------------
optim_func <- function(variable.parameters) {

}

out <- optim(
  par = c(0.025, 0.025),
  fn = get.rel.cost,
  method = "L-BFGS-B",
  lower = variable.parameters$lower_bounds,
  upper = variable.parameters$upper_bounds,
  control = list(fnscale = -1)
) # turn into a maximization problem
