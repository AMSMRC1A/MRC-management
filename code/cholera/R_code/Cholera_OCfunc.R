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

# replicate oc analysis across multiple parameters
# control_type: character to indicate the type of control to be implemented,
#                   "unique": optimize control uniquely in each patch
#                   "uniform": optimize s.t. control is equal in both patches
#                   "max": keep control at maximum value for entire period
#                   "none": no control for entire period
# return_type: vector to indicate what to return,
#                   "v": return time series of each vaccination
#                   "j": return J values broken down by cases/vacc in each patch
#                   "X": return states


# Apply optimal control---------------------------------------------------------
apply_oc <- function(change_params, guess_v1, guess_v2, 
                     guess_u1, guess_u2, init_x, 
                     ode_fn, bounds, adj_fn, control_type,
                     times, params, tol, return_type) {
  # update parameters
  new_params <- params
  p_loc <- match(names(change_params), names(new_params))
  new_params[p_loc[!is.na(p_loc)]] <- change_params[!is.na(p_loc)]
  if (control_type %in% c("unique", "uniform")) {
    out <- run_oc(
      guess_v1, guess_v2, guess_u1, guess_u2, init_x, bounds, ode_fn, adj_fn,
      times, new_params, tol, control_type
    )
  } else if (control_type %in% c("max", "none")) {
    out <- run_no_optim(bounds, init_x, times, ode_fn, new_params, control_type)
  }
  # for now, return v1, v2 time series and j (in list form)
  ret <- list()
  # KD: could these "ifs" be vectorized?
  if ("v" %in% return_type) {
    # return time series of each vaccination strategy
    ret[["v"]] <- cbind(time = times, v1 = out$v1, v2 = out$v2)
  }
  if ("u" %in% return_type){
    #return time series of eahc sanitation strategy
    ret[["u"]] <- cbind(time = times, u1 = out$u1, u2 = out$u2)
  }
  if ("j" %in% return_type) {
    # return j values broken down by cases/vacc in each patch
    ret[["j"]] <- out$j
  }
  if ("X" %in% return_type) {
    # return states
    ret[["X"]] <- out$x
  }
  return(list(ret))
}


# function to run ode and calculate j without optimizing
run_no_optim <- function(bounds, init_x, times, ode_fn, params, control_type) {
  if (control_type == "none") {
    params$v1 <- 0
    params$v2 <- 0
    params$u1 <- 0
    params$u2 <- 0
  } 
  else if (control_type == "max") {
    params$v1 <- bounds["M1_max"]
    params$v2 <- bounds["M2_max"]
    params$u1 <- bounds["U1_max"]
    params$u2 <- bounds["U2_max"]
  }
  out <- ode(y = init_x, times = times, func = ode_fn, parms = params)
  out <- as.data.frame(out)
  j <- calc_j(times, out, params)
  return(list(x = out, v1 = params$v1, v2 = params$v2, u1 = params$u1, u2 = params$u2, j = j))
}

# function to implement optimal control analysis
run_oc <- function(guess_v1, guess_v2, 
                   guess_u1, guess_u2,
                   init_x, bounds, ode_fn, adj_fn,
                   times, params, tol, control_type) {
  # setup variables
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
  return(oc)
}

# implement loop for optimization
oc_optim <- function(v1, v2, u1, u2, x, lambda, # initial guesses
                     IC, lambda_init, # state ICs & final time adjoints
                     bounds, tol, ode_fn, adj_fn,
                     times, params, control_type) {
  with(as.list(bounds), {
    counter <- 1
    counter_limit <- 50
    test <- -1
    while (test < 0 & counter < counter_limit) {
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
      # calculate v1* and v2*
      temp_v <- calc_opt_v(params, lambda, x, control_type)
      temp_u <- calc_opt_u(params, lambda, x, control_type)
      # include bounds
      v1 <- pmin(M1_max, pmax(M1_min, temp_v$temp_v1))
      v2 <- pmin(M2_max, pmax(M2_min, temp_v$temp_v2))
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
      print(paste(counter,test))
    }
    if (counter == counter_limit) {warning("Hit the limit without convergence!", immediate. = TRUE)}
    return(list(
      x = x, lambda = lambda, v1 = v1, v2 = v2, u1 = u1, u2 = u2,
      j = calc_j(times, cbind(as.data.frame(x), v1 = v1, v2 = v2, u1 = u1, u2 = u2), params)
    ))
  })
}

calc_opt_v <- function(params, lambda, x, control_type) {
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

calc_opt_u <- function(params, lambda, x, control_type){
  if (control_type == "uniform"){
    temp_u1 <- ((lambda[,"lambda2"] + params$b1 - lambda[,"lambda1"]) * params$beta_W1 * x[,"S1"] * x[,"W1"] - params$D1 + 
      (lambda[,"lambda6"] + params$b2 - lambda[,"lambda5"]) * params$beta_W2 * x[,"S2"] * x[,"W2"] - params$D2)/
      (2*(params$eta1 + params$eta2))
    temp_u2 <- temp_u1
  }
  if (control_type == "unique"){
    temp_u1 <- ((params$beta_W1 * x[, "S1"] * x[,"W1"])*(params$b1 - lambda[, "lambda1"] + lambda[, "lambda2"]) - params$D1)/
                  (2 * params$eta1)
    temp_u2 <- ((params$beta_W2 * x[, "S2"] * x[,"W2"])*(params$b2 - lambda[, "lambda5"] + lambda[, "lambda6"]) - params$D2)/
      (2 * params$eta2)
  }
  return(list(temp_u1 = temp_u1, temp_u2 = temp_u2))
}

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


calc_j <- function(times, optim_states, params) {
  x <- times
  j_ints <- list(
    case1 = expression(b1 * (beta_I1 * S1 * I1 + (1-u1) * beta_W1 * S1 * W1)),
    case2 = expression(b2 * (beta_I2 * S2 * I2 + (1-u2) * beta_W2 * S2 * W2)),
    vacc1 = expression(C1 * v1 * S1 + epsilon1 * (v1^2)),
    vacc2 = expression(C2 * v2 * S2 + epsilon2 * (v2^2)),
    sani1 = expression(D1 * u1 + eta1 * u1^2), 
    sani2 = expression(D2 * u2 + eta2 * u2^2)
  )
  j_vals <- lapply(j_ints, function(x) {
    apply(optim_states, 1, eval_j_integrand, params = params, integrand = x)
  })
  return(data.frame(
    j_case1 = trapz(x, j_vals[[1]]),
    j_case2 = trapz(x, j_vals[[2]]),
    j_vacc1 = trapz(x, j_vals[[3]]),
    j_vacc2 = trapz(x, j_vals[[4]]),
    j_sani1 = trapz(x, j_vals[[5]]),
    j_sani2 = trapz(x, j_vals[[6]])
  ))
}

