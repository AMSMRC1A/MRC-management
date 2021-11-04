## R resource
# http://desolve.r-forge.r-project.org/user2014/tutorial.pdf (see "Forcing" section)
# https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf

## Pseudo-code
# initial guesses, conditions
# while
#   solve states
#     interpolate control
#   solve adjoint
#     interpolate states and control over time
#   calculate optimal control
#   test optimal control

# replicate optimal control analysis across multiple parameter sets
# control_type: character to indicate the type of control to be implemented,
#                   "unique": optimize control uniquely in each patch
#                   "uniform": optimize s.t. control is equal in both patches
#                   "max": keep control at maximum value for entire period
#                   "none": no control for entire period
# return_type: vector to indicate what to return,
#                   "v": return time series of each vaccination
#                   "j": return J values broken down by cases/vacc in each patch
#                   "X": return states
#                   "all": return data as in "v", "j", and "x"


# Apply optimal control---------------------------------------------------------

# Helper function 'param_changer':
# Update "params" using only the variable names in "change_params"
param_changer <- function(change_params, params) {
  new_params <- params
  p_loc <- match(names(change_params), names(new_params))
  new_params[p_loc[!is.na(p_loc)]] <- change_params[!is.na(p_loc)]
  ret(new_params)
}

# Function 'apply_oc':
# Change parameters (according to "change_params") then apply optimal control
# KD: update this to remove the 'param_changer' call. That can be done prior
#     to calling this function. Reduces the num. of parameters in this function.
#     Consider renaming this function to make it harder to confuse with "run_oc"
apply_oc <- function(change_params, guess_v1, guess_v2, init_x, bounds,
                     ode_fn, adj_fn, control_type,
                     times, params, tol, return_type) {

  # update parameters
  new_params <- param_changer(change_params, params)

  # update settings to match the control type
  if (control_type %in% c("unique", "uniform")) {
    out <- run_oc(
      guess_v1, guess_v2, init_x, bounds, ode_fn, adj_fn,
      times, new_params, tol, control_type
    )
  } else if (control_type %in% c("max", "none")) {
    out <- run_no_optim(bounds, init_x, times, ode_fn, new_params, control_type)
  }

  # Return according to the chosen return_type
  ret <- list()
  # KD: could these "ifs" be vectorized?
  # return_type "all": return all of the below
  if (return_type == "all") {
    return_type <- c("v", "j", "X")
  } # cheap shortcut
  # return_type "v": return time series of each vaccination strategy
  if ("v" %in% return_type) {
    ret[["v"]] <- cbind(time = times, v1 = out$v1, v2 = out$v2)
  }
  # return_type "j": return j values broken down by cases/vacc in each patch
  if ("j" %in% return_type) {
    ret[["j"]] <- out$j
  }
  # return_type "X": return time series of each compartment
  if ("X" %in% return_type) {
    ret[["X"]] <- out$x
  }

  return(list(ret))
}

# Function 'run_no_optim':
# Calculate time series and cost, j, with no optimization
run_no_optim <- function(bounds, init_x, times, ode_fn, params, control_type) {
  if (control_type == "none") {
    params$v1 <- 0
    params$v2 <- 0
  } else if (control_type == "max") {
    params$v1 <- bounds$M1
    params$v2 <- bounds$M2
  }
  out <- ode(y = init_x, times = times, func = ode_fn, parms = params)
  out <- as.data.frame(out)
  j <- calc_j(times, out, params)
  return(list(x = out, v1 = params$v1, v2 = params$v2, j = j))
}

# Function 'run_oc':
# run optimal control analysis
run_oc <- function(guess_v1, guess_v2, init_x, bounds, ode_fn, adj_fn,
                   times, params, tol, control_type) {
  # initialize variables
  x <- matrix(0, nrow = length(times), ncol = 9)
  lambda <- matrix(0, nrow = length(times), ncol = 9)
  v1 <- guess_v1
  v2 <- guess_v2
  if (control_type == "uniform") {
    v2 <- v1
  }
  # final time adjoints
  lambda_init <- rep(0, 8)
  names(lambda_init) <- paste0("lambda", 1:8)
  # implement optimization
  oc <- oc_optim(
    v1, v2, x, lambda,
    IC, lambda_init,
    bounds, tol, ode_fn, adj_fn,
    times, params, control_type
  )
  return(oc)
}

# Sub-function 'oc_optim' (used in 'run_oc'):
# Used in 'run_oc' for the optimization loop
oc_optim <- function(v1, v2, x, lambda, # initial guesses
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
      oldx <- x
      oldlambda <- lambda
      # define interpolating functions for v
      v1_interp <- approxfun(times, v1, rule = 2)
      v2_interp <- approxfun(times, v2, rule = 2)
      # solve states
      x <- ode(
        y = IC, times = times, func = ode_fn, parms = params,
        v1_interp = v1_interp, v2_interp = v2_interp
      )
      # define interpolating functions for x (states)
      x_interp <- lapply(2:ncol(x), function(i) {
        approxfun(x[, c(1, i)], rule = 2)
      })
      # solve adjoint equations (backwards)
      lambda <- ode(
        y = lambda_init, times = rev(times), func = adj_fn, parms = params,
        v1_interp = v1_interp, v2_interp = v2_interp, x_interp = x_interp, x = x
      )
      lambda <- lambda[nrow(lambda):1, ]
      # calculate v1* and v2*
      temp <- calc_opt_v(params, lambda, x, control_type)
      # include bounds
      v1 <- pmin(M1, pmax(0, temp$temp_v1))
      v2 <- pmin(M2, pmax(0, temp$temp_v2))
      # update control
      v1 <- 0.5 * (v1 + oldv1)
      v2 <- 0.5 * (v2 + oldv2)
      # recalculate test
      test <- min(
        tol * norm_oc(c(v1, v2)) - norm_oc(c(oldv1, oldv2) - c(v1, v2)),
        tol * norm_oc(x[, -1]) - norm_oc(oldx[, -1] - x[, -1]),
        tol * norm_oc(lambda[, -1]) - norm_oc(oldlambda[, -1] - lambda[, -1])
      )
      counter <- counter + 1
    }
    return(list(
      x = x, lambda = lambda, v1 = v1, v2 = v2,
      j = calc_j(times, cbind(as.data.frame(x), v1 = v1, v2 = v2), params)
    ))
  })
}

# Sub-function 'calc_opt_v' (used in 'oc_optim'):
# Calculate optimal vaccination strategy
calc_opt_v <- function(params, lambda, x, control_type) {
  params <- as_tibble(as.list(params))
  # KD: !!! for consistency, change "lambda" to use same entry calls as params
  if (control_type == "uniform") {
    temp_v1 <- (((lambda[, "lambda1"] - params$C1 - lambda[, "lambda3"]) * x[, "S1"]) +
      ((lambda[, "lambda5"] - params$C2 - lambda[, "lambda7"]) * x[, "S2"])) /
      (2 * params$epsilon1 + 2 * params$epsilon2)
    temp_v2 <- temp_v1
  }
  if (control_type == "unique") {
    temp_v1 <- ((lambda[, "lambda1"] - params$C1 - lambda[, "lambda3"]) * x[, "S1"]) /
      (2 * params$epsilon1)
    temp_v2 <- ((lambda[, "lambda5"] - params$C2 - lambda[, "lambda7"]) * x[, "S2"]) /
      (2 * params$epsilon2)
  }
  return(list(temp_v1 = temp_v1, temp_v2 = temp_v2))
}

# Helper function 'norm_oc':
# Re-implments the norm(X,1) command from matlab
norm_oc <- function(x) {
  sum(abs(x))
}

# Helper function 'eval_j_integrand':
# define cost function (j values)
eval_j_integrand <- function(params, optim_states, integrand) {
  with(as.list(c(optim_states, params)), {
    eval(parse(text = integrand))
  })
}

# Sub-function 'calc_j' (used in oc_optim):
# Calculates costs (j) from parameters and oc time series
calc_j <- function(times, optim_states, params) {
  x <- times
  j_ints <- list(
    case1 = expression(b1 * (beta_I1 * S1 * I1 + beta_W1 * S1 * W1)),
    case2 = expression(b2 * (beta_I2 * S2 * I2 + beta_W2 * S2 * W2)),
    vacc1 = expression(C1 * v1 * S1 + epsilon1 * (v1^2)),
    vacc2 = expression(C2 * v2 * S2 + epsilon2 * (v2^2))
  )
  j_vals <- lapply(j_ints, function(x) {
    apply(optim_states, 1, eval_j_integrand, params = params, integrand = x)
  })
  return(data.frame(
    j_case1 = trapz(x, j_vals[[1]]),
    j_case2 = trapz(x, j_vals[[2]]),
    j_vacc1 = trapz(x, j_vals[[3]]),
    j_vacc2 = trapz(x, j_vals[[4]])
  ))
}
