# Function 'run_no_optim':
# Calculate time series and cost, j, with no optimization
run_no_optim <- function(bounds, init_x, times, ode_fn, params, control_type) {
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
}

# Function 'run_oc':
# function to implement optimal control analysis
run_oc <- function(guess_v1, guess_v2, 
                   guess_u1, guess_u2,
                   init_x, bounds, ode_fn, adj_fn,
                   times, params, tol, control_type) {
  # initialize variables
  x <- matrix(0, nrow = length(times), ncol = 9)
  lambda <- matrix(0, nrow = length(times), ncol = 9)
  v1 <- guess_v1
  v2 <- guess_v2
  u1 <- guess_u1
  u2 <- guess_u2
  if (control_type == "uniform") {
    v2 <- v1
    u2 <- u1
  }
  # final time adjoints
  lambda_init <- rep(0, 8)
  names(lambda_init) <- paste0("lambda", 1:8)
  # implement optimization
  oc <- oc_optim(
    v1, v2, u1, u2, x, lambda,
    IC, lambda_init,
    bounds, tol, ode_fn, adj_fn,
    times, params, control_type
  )
  # "oc" is a list with entries:
  #     trajectories = all time-varying values (state variables, controls,
  #                    and adjoints)
  #     j = costs
  return(oc)
}

# Sub-function 'oc_optim' (used in 'run_oc'):
# Used in 'run_oc' for the optimization loop
# implement loop for optimization
oc_optim <- function(v1, v2, u1, u2, x, lambda, # initial guesses
                     IC, lambda_init, # state ICs & final time adjoints
                     bounds, tol, # optimal control settings
                     ode_fn, adj_fn, # ode and adjoint functions
                     times, params, control_type # ode and oc parameters
) {
  with(as.list(bounds), {
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
        y = IC, times = times, func = ode_fn, parms = params,
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
        tol * norm_oc(c(u1, u2)) - norm_oc(c(oldu1, oldu2) - c(u1, u2)), #EH: ASK SUZANNE ABOUT THIS
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

# Utilities --------------------------------------------------------------------
# define norm(X,1) command from matlab
norm_oc <- function(x) {
  sum(abs(x))
}

# define cost function (j values)
eval_j_integrand <- function(params, optim_states, integrand) {
  with(as.list(c(optim_states, params)), {
    eval(parse(text = integrand))
  })
}
