# optimal control --------------------------------------------------------------

#' implement optimal control optimization
#' 
#' @param model string to indicate which model is being run
#' either \code{"uniform"} for the same control being applied in both patches 
#' or \code{"unique"} where control can vary across patches 
#' @param filepath string specifying relative path location of model specific code
#' (i.e., "cholera" or "ebola" folder). Defaults to "models" folder 
#' (assumes MRC-management/code/R) is working directory
#' @param change_params data frame with parameter values to change from baseline
#' @param counter_max integer maximum number of iterations for optimization
#' 
#' #' @return list containing \code{trajectories}, a data.frame with all 
#' time-varying values (state variables, controls, and adjoints), and 
#' \code{j}, a double of the total cost
oc_optim <- function(model, filepath = "models/", 
                     change_params = NA, counter_max = 50) {
  # load baseline objects
  setup <- setup_model(model, filepath)
  with(setup, {
    # update parameters if !is.na(change_params)
    if(!any(is.na(change_params))){
      params <- param_changer(change_params, params)
      IC <- IC_setup(params = params,
                     model_ODE = ode_fn, 
                     IC_init = init_x_uncontrol, 
                     response_time = response_time)
      init_x <- IC$IC
    }
    counter <- 1
    test <- -1
    test_vals <- list() # list to store test
    while (test < 0 & counter < counter_max) {
      # If we hit the max number of iterations, give a chance to increase iterations or quit out
      # if(counter == counter_max - 1){
      #   warning("max number of iterations reached")
      #   new_counter_max <- readline("Input additional number of iterations to perform: ")
      #   if(new_counter_max %in% c("F", "N", "no")) {break} 
      #   else {
      #     new_counter_max <- as.integer(new_counter_max)
      #     counter_max <- counter_max + new_counter_max}
      # }
      # else {
        # set previous control, state, and adjoint
        old_vals <- set_old_variables(c(controls,
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
                                            optimal_control_fn = optimal_control_fn,
                                            old_controls = old_vals)
        # recalculate test
        test <- calc_test_fn(params$tol, controls, x, lambda, old_vals)
        test_vals[[counter]] <- test
        counter <- counter + 1
      }
    # }
    
    trajectories <- cbind(x_df, do.call(cbind, controls))
    j_vals <- calc_j(times, 
                     cbind(as.data.frame(x), 
                           do.call(cbind, controls)),
                     params)
    if(counter == counter_max){
      warning("max number of iterations reached")
      return(list(
        trajectories = NA,
        j = NA,
        uncontrolled = IC$uncontrolled,
        test_vals = unlist(test_vals),
        n_iterations = counter
      ))
    }
    else{
      return(list(
        trajectories = trajectories,
        j = j_vals, 
        uncontrolled = IC$uncontrolled,
        test_vals = unlist(test_vals), 
        n_iterations = counter
      ))
    }
  })
}

#' run optimal control optimization over different parameter sets
#' 
#' @param model string to indicate which model is being run
#' either \code{"cholera"} for the cholera model or \code{"ebola"} for the ebola
#' model
#' @param param_sets data frame composed of parameter values to change from
#'  baseline. Each row corresponds to a different 'experiment' to run
#' 
#' #' @return list containing \code{trajectories}, a data.frame with all 
#' time-varying values (state variables and controls) indexed by the row number
#' from param_sets, and \code{j}, a vector of the costs index by the row number
#' from param_sets

mult_oc_optim <- function(model, param_sets = NA) {
  # Set up loop over param_sets
  # Run loop over param_sets
  vary_params <- foreach(
    i = 1:nrow(param_sets),
    .packages = c("deSolve", "tidyverse", "pracma"),
    # explicitly give 'foreach' the functions and data it needs
    .export = c("param_changer","oc_optim", "setup_model", "set_old_variables",
                "define_interp_fns", "update_optimal_solution")
  ) %dopar% {
    
    oc_optim(model, param_sets[i, ])
  }
  
  # Organize output from loops 
  param_sets$test_case <- 1:nrow(param_sets)
  trajectories <- lapply(1:nrow(param_sets), function(i) {
    return(data.frame(test_case = i, vary_params[[i]]$trajectories))
  }) 
  trajectories <- as_tibble(do.call(rbind, trajectories))
  out <- left_join(param_sets, trajectories)
  j_vals <- lapply(1:nrow(param_sets), function(i) {
    return(data.frame(test_case = i, vary_params[[i]][["j"]]))
  })
  j_vals <- as_tibble(do.call(rbind, j_vals))
  j_vals$j <- apply(j_vals[, c("j_case1", "j_case2", "j_vacc1", "j_vacc2")], 1, sum)
  out <- left_join(out, j_vals)
  return(out)  
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

#' input all model-specific variables for optimal control analysis
#' 
#' @param model string to indicate which values to indicate which model to input
#' @param filepath string specifying relative path location of model specific code
#' (i.e., "cholera" or "ebola" folder). 
#' Defaults to "models" folder (assumes MRC-management/code/R) is working directory
#' 
#' @return list of initial guesses, initial conditions, optimal control and model
#' settings, and ode and adjoint functions
setup_model <- function(model, filepath = "models/"){
  # source files where model-specific information is stored
  source(file.path(paste0(filepath,model,"/",model,"_baseline_params.R")))
  source(file.path(paste0(filepath,model,"/",model,"_functions.R")))
  # define "dictionary" for model setup
  n_states = length(get(paste0("IC_",model,"_uncontrol")))
  # final time adjoints
  lambda_init <- rep(0, n_states)
  names(lambda_init) <- paste0("lambda", 1:n_states)
  # create list of all model-specific objects to load
  setup <- list(
    # initial guesses
    controls = list(
      v1 = guess_v1, 
      v2 = guess_v2, 
      u1 = guess_u1,
      u2 = guess_u2), 
    x = matrix(0, nrow = length(get(paste0("times_",model))), ncol = n_states+1),
    lambda = matrix(0, nrow = length(get(paste0("times_",model))), ncol = n_states+1),
    # ICs for uncontrolled ode
    init_x_uncontrol = get(paste0("IC_",model,"_uncontrol")),
    response_time = get(paste0("response_time_", model)),
    lambda_init = lambda_init, 
    # functions
    ode_fn = get(paste0("ode_",model)), 
    adj_fn = get(paste0("adjoint_",model)),
    optimal_control_fn = get(paste0("optimal_controls_",model)),
    calc_test_fn = get(paste0("calc_test_",model)),
    calc_j = get(paste0("calc_j_",model)),
    # model settings
    times = get(paste0("times_",model)), 
    params =get(paste0("params_",model))
  )
  return(setup)
}

## EH: THERE MAY BE A SIMPLER WAY TO IMPLEMENT THIS
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

#' calculate optimal solution
#' 
#' calculate optimal solution using optimal control characterization, and 
#' update to account for bounds
#' 
#' @param params vector of model and optimal control parameters
#' @param lambda matrix of adjoints over time
#' @param x matrix of states over time
#' either \code{"uniform"} for the same control being applied in both patches 
#' or \code{"unique"} where control can vary across patches 
#' @param optimal_control_fn function to calculate optimal control 
#' characterization
#' @param old_controls list of controls from previous iteration
update_optimal_solution <- function(params, lambda, x,
                                    optimal_control_fn, old_controls){
  with(params, {
    # calculate v1*, v2*, u1*, u2* (or other optimal controls)
    temp_controls <- optimal_control_fn(params, lambda, x, control_type)
    controls <- list()
    for(i in names(temp_controls)){
      # include bounds
      controls[[i]] <- pmin(eval(as.name(paste0(i,"_max"))), 
                            pmax(eval(as.name(paste0(i,"_min"))), 
                                 temp_controls[[i]]))
      # update control
      controls[[i]] <- 0.5 * (controls[[i]] + old_controls[[paste0("old",i)]])
    }
    return(controls)
  })
}

#' Helper function 'param_changer'
#' 
#' Update the \code{params} data frame with the parameter values listed in the
#' \code{change_params} data frame
#' 
#' @param change_params data frame containing new parameter values
#' @param params data frame containing original parameter values
#' 
#' @return new_params data frame with updated parameter values
param_changer <- function(change_params, params) {
  new_params <- params
  p_loc <- match(names(change_params), names(new_params))
  new_params[p_loc[!is.na(p_loc)]] <- change_params[!is.na(p_loc)]
  return(new_params)
}

#' Helper function 'IC_setup'
#' 
#' Determine initial conditions for optimal control by running model 
#' uncontrolled for a given number of days 
#' 
#' @param params vector of parameters to pass into ODE
#' @param model ODE function for given model
#' @param IC_init initial conditions for uncontrolled ODE
#' @param stepsize double indicating time increments for ODE
#' @param response_time time after which to start optimal control
#' 
#' @return vector of optimal control initial conditions
IC_setup <- function(params, model_ODE, IC_init, stepsize = 0.1, response_time){
  uncontrolled <- ode(
    y = IC_init,
    times = seq(1,response_time, stepsize),
    func = model_ODE,
    parms = params
  )
  # set IC based on response_time
  IC <- as.double(uncontrolled[uncontrolled[, "time"] == response_time, -1])
  names(IC) <- colnames(uncontrolled[, -1])
  return(list(uncontrolled = uncontrolled, IC = IC))
}