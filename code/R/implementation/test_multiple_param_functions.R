# Helper function 'param_changer':
# Update the "params" data frame with the parameter values listed in the
# "change_params" data frame
param_changer <- function(change_params, params) {
  new_params <- params
  p_loc <- match(names(change_params), names(new_params))
  new_params[p_loc[!is.na(p_loc)]] <- change_params[!is.na(p_loc)]
  return(new_params)
}

# Function 'apply_oc':
# Change parameters (according to "change_params") then apply optimal control
# KD to do: update this to remove the 'param_changer' call. That can be done
#           prior to calling this function. Reduces the num. of parameters in
#           this function. Consider renaming this function to make it harder to
#           confuse with "run_oc"
apply_oc <- function(change_params, guess_v1, guess_v2,
                     guess_u1, guess_u2, init_x, bounds,
                     ode_fn, adj_fn, control_type,
                     times, params, tol) {
  
  # update parameters
  new_params <- param_changer(change_params, params)
  
  # update settings to match the control type input
  if (control_type %in% c("unique", "uniform")) {
    out <- run_oc(
      guess_v1, guess_v2, guess_u1, guess_u2, init_x, bounds, ode_fn, adj_fn,
      times, new_params, tol, control_type
    )
  } else if (control_type %in% c("max", "none")) {
    out <- run_no_optim(bounds, init_x, times, ode_fn, new_params, control_type)
  }
  # got rid of return_types:
  # "v" data is now included in x, calculating j values not computationally expensive
  return(out)
}

#### FUNCTIONS TO IMPLEMENT OPTIMAL CONTROL ANALYSES ####

# run OC over multiple parameter sets-------------------------------------------
# test_params: data.frame containing the parameters and values to vary
#              including control_type, each row represents one case to be tested
#
# base_params: vector containing baseline params to be used if not varying
#              includes both biological and OC params
test_mult_params <- function(test_params, base_params,
                             guess_v1, guess_v2, 
                             guess_u1, guess_u2, 
                             IC, bounds, times, tol) {
  # Run optimal control calculations across test_params dataframe
  vary_params <- foreach(
    i = 1:nrow(test_params),
    .packages = c("deSolve", "tidyverse", "pracma"),
    # explicitly give 'foreach' the functions and data it needs
    .export = c("apply_oc","param_changer","run_oc","oc_optim", "run_no_optim",
                "adj","calc_opt_v", "calc_opt_u",
                "norm_oc","calc_j","eval_j_integrand",
                "chol","params", "oc_params")
  ) %dopar% {
    apply_oc(
      change_params = test_params[i, ],
      guess_v1 = guess_v1, guess_v2 = guess_v2,
      guess_u1 = guess_u1, guess_u2 = guess_u2,
      init_x = IC, bounds = bounds,
      ode_fn = chol, adj_fn = adj,
      control_type = test_params[i, "control_type"],
      times = times, params = base_params,
      tol = tol
    )
  }
  test_params$test_case <- 1:nrow(test_params)
  states <- lapply(1:nrow(test_params), function(i) {
    return(data.frame(test_case = i, vary_params[[i]]$trajectories))
  })
  states <- as.data.frame(do.call(rbind, states))
  states <- left_join(test_params, states)
  j_vals <- lapply(1:nrow(test_params), function(i) {
    return(data.frame(test_case = i, vary_params[[i]][["j"]]))
  })
  j_vals <- as.data.frame(do.call(rbind, j_vals))
  j_vals <- left_join(test_params, j_vals)
  j_vals$j <- apply(j_vals[, c("j_case1", "j_case2", "j_vacc1", "j_vacc2")], 1, sum)
  return(list(states = states, j_vals = j_vals))
}