#### FUNCTIONS TO IMPLEMENT OMPTIMAL CONTROL ANALYSES ####

# run OC over multiple parameter sets-------------------------------------------
# test_params: data.frame containing the parameters and values to vary
#              including control_type, each row represents one case to be tested
# return_type: vector containing "v", "j", "X" (see apply_oc() in OCfunc.R for
#              full definitions)
# base_params: vector containing baseline params to be used if not varying
#              includes both biological and OC params
test_mult_params <- function(test_params, return_type, base_params,
                             guess_v1, guess_v2, 
                             guess_u1, guess_u2, 
                             IC, bounds, times, tol) {
  # Run optimal control calculations across test_params dataframe
  vary_params <- foreach(
    i = 1:nrow(test_params),
    .packages = c("deSolve", "tidyverse", "pracma"),
    # explicitly give 'foreach' the functions and data it needs
    .export = c(
      "apply_oc", "guess_v1", "guess_v2",
      "guess_u1", "guess_u2", "IC",
      "bounds", "chol", "adj", "times",
      "params", "oc_params", "run_oc", "tol",
      "oc_optim", "calc_opt_v", "calc_opt_u", "norm_oc", "calc_j",
      "eval_j_integrand"
    ),
    .combine = rbind
  ) %dopar% {
    apply_oc(
      change_params = test_params[i, ],
      guess_v1 = guess_v1, guess_v2 = guess_v2,
      guess_u1 = guess_u1, guess_u2 = guess_u2,
      init_x = IC, bounds = bounds,
      ode_fn = chol, adj_fn = adj,
      times = times, params = c(params, oc_params),
      tol = tol,
      control_type = test_params[i, "control_type"],
      return_type = return_type
    )
  }
  # vary_params <- do.call(rbind, vary_params)
  return(vary_params)
}

# reformat list output to data.frames-------------------------------------------
reformat_mult_params_output <- function(output, test_params) {
  # Assign row numbers as a column "test_case"
  test_params$test_case <- 1:nrow(test_params)
  # Invoke individual reformatting functions for each output type
  # use first element in output list to determine types of output included
  reformatted <- lapply(
    names(output[[1]]),
    function(i) {
      return(get(paste0("reformat_output_", i))(output, test_params))
    }
  )
  names(reformatted) <- names(output[[1]])
  return(reformatted)
}

reformat_output_X <- function(output, test_params) {
  # reformat state vector time series
  states <- lapply(1:length(output), function(i) {
    return(data.frame(test_case = i, output[[i]][["X"]]))
  })
  states <- as.data.frame(do.call(rbind, states))
  states <- left_join(test_params, states)
  return(states)
}

reformat_output_v <- function(output, test_params) {
  # reformat optimal control strategy time series
  v_timeseries <- lapply(1:length(output), function(i) {
    return(data.frame(test_case = i, output[[i]][["v"]]))
  })
  v_timeseries <- as.data.frame(do.call(rbind, v_timeseries))
  v_timeseries <- left_join(test_params, v_timeseries)
  return(v_timeseries)
}

reformat_output_u <- function(output, test_params) {
  # reformat optimal control strategy time series
  u_timeseries <- lapply(1:length(output), function(i) {
    return(data.frame(test_case = i, output[[i]][["u"]]))
  })
  u_timeseries <- as.data.frame(do.call(rbind, u_timeseries))
  u_timeseries <- left_join(test_params, u_timeseries)
  return(u_timeseries)
}

reformat_output_j <- function(output, test_params) {
  # reformat cost calculations
  j_vals <- lapply(1:length(output), function(i) {
    return(data.frame(test_case = i, output[[i]][["j"]]))
  })
  j_vals <- as.data.frame(do.call(rbind, j_vals))
  j_vals <- left_join(test_params, j_vals)
  j_vals$j <- apply(j_vals[, c("j_case1", "j_case2", "j_vacc1", "j_vacc2")], 1, sum)
  return(j_vals)
}
