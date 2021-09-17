## R resource
# http://desolve.r-forge.r-project.org/user2014/tutorial.pdf (see "Forcing" section)
# https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf

## Psuedo-code
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
apply_oc = function(change_params,guess_v1, guess_v2, init_x, bounds,
                    ode_fn, adj_fn, control_type,
                    times, params, delta) {
  # update parameters
  new_params <- params 
  p_loc <- match(names(change_params), names(new_params))
  new_params[p_loc[!is.na(p_loc)]] = change_params[!is.na(p_loc)]
  if(control_type %in% c("unique", "uniform")){
    out <- run_oc(guess_v1, guess_v2, init_x, bounds, ode_fn, adj_fn,
                  times, new_params, delta, control_type)
  }
  else if(control_type %in% c("max", "none")){
    out <- run_no_optim(bounds, init_x, times, ode_fn, new_params, control_type)
  }
  # for now, return v1, v2 time series and j (in list form)
  ret <- list(list(ts = cbind(time = times, v1 = out$v1, v2 = out$v2), j = out$j))
  return(ret)
}


# function to run ode and calculate j without optimizing
run_no_optim = function(bounds, init_x, times, ode_fn, params, control_type){
  if(control_type == "none"){
    params$v1 = 0
    params$v2 = 0
  }
  else if(control_type == "max"){
    params$v1 = bounds[1]
    params$v2 = bounds[2]
  }
  out <- ode(y = init_x, times = times, func = ode_fn, parms = params)
  out = as.data.frame(out)
  j <- calc_j(times,out, params)
  return(list(x = out, v1 = params$v1, v2 = params$v2, j = j))
}

# function to implement optimal control analysis
run_oc = function(guess_v1, guess_v2, init_x, bounds,ode_fn, adj_fn,
                  times, params, delta, control_type){
  # setup variables 
  x = matrix(0, nrow = length(times), ncol = 9)
  lambda = matrix(0, nrow = length(times), ncol = 9)
  v1 = guess_v1
  v2 = guess_v2
  if(control_type == "uniform"){
    v2 = v1
  }
  # final time adjoints
  lambda_init = rep(0,8);
  names(lambda_init) = paste0("lambda",1:8)
  # implement optimization
  oc = oc_optim(v1, v2, x, lambda, 
                IC, lambda_init, 
                bounds, delta, ode_fn, adj_fn, 
                times, params, control_type)
  oc$j <- calc_j(times,cbind(oc$x, v1 = oc$v1, v2 = oc$v2), params)
  # oc$j <- calc_j(params = params, 
  #             optim_states = cbind(oc$x, v1 = oc$v1, v2 = oc$v2), 
  #             integrand_fn = j_integrand, 
  #             lower_lim = min(times), upper_lim = max(times), step_size = (range(times)[2] - range(times)[1])/(length(times)-1))
  return(oc)
}

# implement loop for optimization
oc_optim = function(v1, v2, x, lambda, # initial guesses
                    IC, lambda_init, # state ICs & final time adjoints
                    bounds, delta, ode_fn, adj_fn,
                    times, params, control_type){
  with(as.list(bounds),{
    counter = 1
    test = -1
    while(test < 0 & counter < 50){
      # set previous control, state, and adjoint 
      oldv1 = v1
      oldv2 = v2
      oldx = x
      oldlambda = lambda
      # define interpolating functions for v
      v1_interp <- approxfun(times, v1, rule = 2)
      v2_interp <- approxfun(times, v2, rule = 2)
      # solve states
      x <- ode(y = IC, times = times, func = ode_fn, parms = params, 
               v1_interp = v1_interp, v2_interp = v2_interp)
      # define interpolating functions for x (states)
      x_interp <- lapply(2:ncol(x), function(i){approxfun(x[,c(1,i)], rule = 2)})
      # solve adjoint equations (backwards)
      lambda <- ode(y = lambda_init, times = rev(times), func = adj_fn, parms = params,
                    v1_interp = v1_interp, v2_interp = v2_interp, x_interp = x_interp, x = x)
      lambda <- lambda[nrow(lambda):1,]
      # calculate v1* and v2*
      temp <- calc_opt_v(params, lambda, x, control_type)
      # include bounds
      v1 = pmin(M1, pmax(0, temp$temp_v1))
      v2 = pmin(M2, pmax(0, temp$temp_v2))
      # update control
      v1 = 0.5*(v1 + oldv1)
      v2 = 0.5*(v2 + oldv2)
      # recalculate test
      test <- min(delta*norm_oc(c(v1,v2))-norm_oc(c(oldv1,oldv2)-c(v1,v2)),
                  delta*norm_oc(x[,-1])-norm_oc(oldx[,-1]-x[,-1]),
                  delta*norm_oc(lambda[,-1])-norm_oc(oldlambda[,-1]-lambda[,-1]))
      print(counter)
      print(test)
      calc_j(times,cbind(as.data.frame(x), v1 = v1, v2 = v2), params)
      counter <- counter + 1
    }
    return(list(x = x, lambda = lambda, v1 = v1, v2 = v2))
  })
}

calc_opt_v <- function(params, lambda, x, control_type){
  if(control_type == "uniform"){
    temp_v1 <- (((lambda[,"lambda1"] - params$C1 - lambda[,"lambda3"])*x[,"S1"]) + 
      ((lambda[,"lambda5"] - params$C2 - lambda[,"lambda7"])*x[,"S2"]))/
      (2*params$epsilon1 + 2*params$epsilon2)
    temp_v2 = temp_v1
  }
  if(control_type == "unique"){
    temp_v1 <- ((lambda[,"lambda1"] - params$C1 - lambda[,"lambda3"])*x[,"S1"])/
      (2*params$epsilon1)
    temp_v2 <- ((lambda[,"lambda5"] - params$C2 - lambda[,"lambda7"])*x[,"S2"])/
      (2*params$epsilon2)
  }
  return(list(temp_v1 = temp_v1, temp_v2 = temp_v2))
}

# define norm(X,1) command from matlab
norm_oc <- function(x){sum(abs(x))}

# define cost function (j values)
eval_j_integrand <- function(params, optim_states, integrand){
  with(as.list(c(optim_states, params)),{
    eval(parse(text = integrand))
  })
}


calc_j <- function(times,optim_states,params){
  x <- times
  j_ints <- list(case1 = expression(b1 * (beta_I1*S1*I1 + beta_W1*S1*W1)),
                 case2 = expression(b2 * (beta_I2*S2*I2 + beta_W2*S2*W2)),
                 vacc1 = expression(C1*v1*S1 + epsilon1*(v1^2)),
                 vacc2 = expression(C2*v2*S2 + epsilon2*(v2^2)))
  j_vals <- lapply(j_ints, function(x){apply(optim_states, 1, eval_j_integrand, params = params, integrand = x)})
  return(data.frame(j_case1 = trapz(x,j_vals[[1]]),
                    j_case2 = trapz(x,j_vals[[2]]),
                    j_vacc1 = trapz(x,j_vals[[3]]),
                    j_vacc2 = trapz(x,j_vals[[4]])))
}


