#### FUNCTIONS TO IMPLEMENT OPTIMAL CONTROL ANALYSES ####

# run OC over multiple parameter sets-------------------------------------------
# test_params: data.frame containing the parameters and values to vary
#              including control_type, each row represents one case to be tested
#
# base_params: vector containing baseline params to be used if not varying
#              includes both biological and OC params
test_mult_params <- function(test_params, base_params,
                             guess_v1, guess_v2, IC, bounds, times, tol) {
  # Run optimal control calculations across test_params dataframe
  vary_params <- foreach(
    i = 1:nrow(test_params),
    .packages = c("deSolve", "tidyverse", "pracma"),
    # explicitly give 'foreach' the functions and data it needs
    .export = c("apply_oc","param_changer","run_oc","oc_optim", "run_no_optim",
                "adj","calc_opt_v","norm_oc","calc_j","eval_j_integrand",
                "chol","params", "oc_params")
  ) %dopar% {
    apply_oc(
      change_params = test_params[i, ],
      guess_v1 = guess_v1, guess_v2 = guess_v2,
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
