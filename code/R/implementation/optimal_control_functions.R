# optimal control --------------------------------------------------------------

#' implement optimal control optimization
#' 
#' @param model string to indicate which model is being run
#' @param control_type character to define the type of control being implemented;
#' either \code{"uniform"} for the same control being applied in both patches 
#' or \code{"unique"} where control can vary across patches 
#' 
#' #' @return list containing \code{trajectories}, a data.frame with all 
#' time-varying values (state variables, controls, and adjoints), and 
#' \code{j}, a double of the total cost
oc_optim <- function(model,  control_type) {
  setup <- setup_model(model)
  with(setup, {
    counter <- 1
    test <- -1
    while (test < 0 & counter < 50) {
      # set previous control, state, and adjoint
      oldv1 <- v1
      oldv2 <- v2
      oldu1 <- u1
      oldu2 <- u2
      oldx <- x
      oldlambda <- lambda
      # define interpolating functions for v
      v1_interp <- approxfun(times, v1, rule = 2)
      v2_interp <- approxfun(times, v2, rule = 2)
      u1_interp <- approxfun(times, u1, rule = 2)
      u2_interp <- approxfun(times, u2, rule = 2)
      # solve states
      x <- ode(
        y = init_x, times = times, func = ode_fn, parms = params,
        v1_interp = v1_interp, v2_interp = v2_interp,
        u1_interp = u1_interp, u2_interp = u2_interp
      )
      x_df <- as.data.frame(x)
      # define interpolating functions for x (states)
      x_interp <- lapply(2:ncol(x), function(i) {
        approxfun(x[, c(1, i)], rule = 2)
      })
      # solve adjoint equations (backwards)
      lambda <- ode(
        y = lambda_init, times = rev(times), func = adj_fn, parms = params,
        v1_interp = v1_interp, v2_interp = v2_interp, 
        u1_interp = u1_interp, u2_interp = u2_interp,
        x_interp = x_interp, x = x
      )
      lambda <- lambda[nrow(lambda):1, ]
      lambda_df <- as.data.frame(lambda)
      # calculate v1*, v2*, u1*, u2*
      temp_v <- calc_opt_v(params, lambda_df, x_df, control_type)
      temp_u <- calc_opt_u(params, lambda_df, x_df, control_type)
      # include bounds
      v1 <- pmin(V1_max, pmax(V1_min, temp_v$temp_v1))
      v2 <- pmin(V2_max, pmax(V2_min, temp_v$temp_v2))
      u1 <- pmin(U1_max, pmax(U1_min, temp_u$temp_u1))
      u2 <- pmin(U2_max, pmax(U2_min, temp_u$temp_u2))
      # update control
      v1 <- 0.5 * (v1 + oldv1)
      v2 <- 0.5 * (v2 + oldv2)
      u1 <- 0.5 * (u1 + oldu1)
      u2 <- 0.5 * (u2 + oldu2)
      # recalculate test
      test <- min(
        tol * norm_oc(c(v1, v2)) - norm_oc(c(oldv1, oldv2) - c(v1, v2)),
        tol * norm_oc(c(u1, u2)) - norm_oc(c(oldu1, oldu2) - c(u1, u2)),
        tol * norm_oc(x[, -1]) - norm_oc(oldx[, -1] - x[, -1]),
        tol * norm_oc(lambda[, -1]) - norm_oc(oldlambda[, -1] - lambda[, -1])
      )
      counter <- counter + 1
    }
    trajectories <- x_df
    trajectories$v1 <- v1
    trajectories$v2 <- v2
    trajectories$u1 <- u1
    trajectories$u2 <- u2
    j_vals <- calc_j(times, 
                     cbind(as.data.frame(x), v1 = v1, v2 = v2, u1 = u1, u2 = u2),
                     params)
    return(list(
      trajectories = trajectories,
      j = j_vals
    ))
  })
}

# no control -------------------------------------------------------------------

#' calculate time series and cost, j, with no optimization
#' 
#' @inheritParams oc_optim
#' 
#' @return list containing \code{trajectories}, a data.frame with all 
#' time-varying values (state variables, controls, and adjoints), and 
#' \code{j}, a double of the total cost
run_no_optim <- function(model, control_type) {
  setup <- setup_model(model)
  with(setup, {
  if (control_type == "none") {
    params$v1 <- 0
    params$v2 <- 0
    params$u1 <- 0
    params$u2 <- 0
  } 
  else if (control_type == "max") {
    params$v1 <- bounds$V1_max
    params$v2 <- bounds$V2_max
    params$u1 <- bounds$U1_max
    params$u2 <- bounds$U2_max
  }
  out <- ode(y = init_x, times = times, func = ode_fn, parms = params)
  trajectories <- as.data.frame(out)
  trajectories$v1 <- params$v1
  trajectories$v2 <- params$v2
  trajectories$u1 <- params$u1
  trajectories$u2 <- params$u2
  j <- calc_j(times, out, params)
  return(list(trajectories = trajectories, j = j))
  })
}


# utilities --------------------------------------------------------------------
#' define norm(X,1) command from matlab
#' 
#' @param x vector
norm_oc <- function(x) {
  sum(abs(x))
}

#' input all model-specific variables for optimal control analysis
#' 
#' @param model string to indicate which values to indicate which model to input
#' 
#' @return list of initial guesses, intial conditions, optimal control and model
#' settings, and ode and adjoint functions
setup_model <- function(model){
  # final time adjoints
  lambda_init <- rep(0, 8)
  names(lambda_init) <- paste0("lambda", 1:8)
  # "dictionary" for cholera setup
  if(model == "cholera"){
    source("models/cholera/cholera_baseline_params.R")
    source("models/cholera/cholera_functions.R")
    # initial guesses
    setup <- list(v1 = guess_v1, 
                  v2 = guess_v2, 
                  u1 = guess_u1, 
                  u2 = guess_u2, 
                  x = matrix(0, nrow = length(times), ncol = 9),
                  lambda = matrix(0, nrow = length(times), ncol = 9),
                  # ICs for ode solver
                  init_x = IC_cholera, 
                  lambda_init = lambda_init, 
                  # optimal control settings
                  tol = tol_cholera,
                  bounds = bounds_cholera, 
                  # functions
                  ode_fn = ode_cholera, 
                  adj_fn = adjoint_cholera, 
                  # model settings
                  times = times_cholera, 
                  params = params_cholera, 
    )
  }
  return(setup)
}
