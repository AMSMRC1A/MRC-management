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
      old_controls <- set_old_variables(c(controls,
                          list(x = x,lambda = lambda)))
      # define interpolating functions for v
      interp_controls <- define_interp_fns(controls, times)
      # solve states
      x <- ode(
        y = init_x, times = times, func = ode_fn, parms = params,
        interp_controls = interp_controls
      )
      x_df <- as.data.frame(x)
      # define interpolating functions for x (states)
      x_interp <- lapply(2:ncol(x), function(i) {
        approxfun(x[, c(1, i)], rule = 2)
      })
      # solve adjoint equations (backwards)
      lambda <- ode(
        y = lambda_init, times = rev(times), func = adj_fn, parms = params,
        interp_controls = interp_controls,
        x_interp = x_interp, x = x
      )
      lambda <- lambda[nrow(lambda):1, ]
      lambda_df <- as.data.frame(lambda)
      # calculate new controls 
      controls <- update_optimal_solution(params = params, 
                                          lambda = lambda_df, 
                                          x = x_df, 
                                          control_type = control_type, 
                                          optimal_control_fn = optimal_control_fn,
                                          bounds = bounds, 
                                          old_controls = old_controls)
      # recalculate test
      test <- calc_test_fn(tol, controls, old_controls, x, lambda)
      counter <- counter + 1
    }
    trajectories <- cbind(x_df, do.call(cbind, controls))
    j_vals <- calc_j(times, 
                     cbind(as.data.frame(x), 
                           do.call(cbind, controls)),
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
    # create list of all model-specific objects to load
    setup <- list(
      # initial guesses
      controls = list(
        v1 = guess_v1, 
        v2 = guess_v2, 
        u1 = guess_u1, 
        u2 = guess_u2), 
      x = matrix(0, nrow = length(times_cholera), ncol = 9),
      lambda = matrix(0, nrow = length(times_cholera), ncol = 9),
      # ICs for ode solver #EH: some of these variable names are confusing
      init_x = IC_cholera, 
      lambda_init = lambda_init, 
      # optimal control settings
      tol = tol_cholera,
      bounds = bounds_cholera, 
      # functions
      ode_fn = ode_cholera, 
      adj_fn = adjoint_cholera,
      optimal_control_fn = optimal_controls_cholera,
      calc_test_fn = calc_test_cholera,
      calc_j = calc_j_cholera,
      # model settings
      times = times_cholera, 
      params = params_cholera
    )
  }
  return(setup)
}

## EH: THERE MAY BE A SIMPLER WAY TO IMPLEMENT THESE

#' add renamed variables to the parent environment
#' 
#' all variables in \code{vars} will generate a new, identical object in the  
#' parent environment as "old" + variable name (e.g., \code{v1} will 
#' generate \code{oldv1}) 
#' 
#' @param vars list of variables to be renamed
set_old_variables <- function(vars){
  renamed_old <- list()
  for(i in names(vars)){
    renamed_old[[paste0("old",i)]] <- vars[[i]]
  }
  return(renamed_old)
}

#' create interpolation functions for controls
#' 
#' all variables in \code{vars} will generate a new object 
#' that is the function to linearly interpolate the values in \code{vars};
#' the name of the new object is variable name + "_interp"
#' (e.g., \code{v1} will generate \code{v1_interp}) 
#' 
#' @param vars list of variables to be renamed
#' @param times vector of time points to interpolate over
#' 
#' @return list of interpolation functions for each element in \code{vars}
define_interp_fns <- function(vars, times){
  interp_fns <- list()
  for(i in names(vars)){
    interp_fns[[i]] <-  approxfun(times, vars[[i]], rule = 2) 
  }
  return(interp_fns)
}

### ADD DOCUMENTATION
update_optimal_solution <- function(params, lambda, x, control_type, 
                                    optimal_control_fn, bounds, old_controls){
  # calculate v1*, v2*, u1*, u2* (or other optimal controls)
  temp_controls <- optimal_control_fn(params, lambda, x, control_type)
  controls <- list()
  for(i in names(temp_controls)){
    # include bounds
    controls[[i]] <- pmin(bounds[[paste0(i,"_max")]], 
                          pmax(bounds[[paste0(i,"_min")]], temp_controls[[i]]))
    # update control
    controls[[i]] <- 0.5 * (controls[[i]] + old_controls[[paste0("old",i)]])
  }
  return(controls)
}

