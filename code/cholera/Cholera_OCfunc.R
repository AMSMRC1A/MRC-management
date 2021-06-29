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

# function to implement optimal control analysis
run_oc = function(guess_v1, guess_v2, init_x, bounds,ode_fn, adj_fn,
                  t, params, oc_params, delta){
  # setup variables 
  x = matrix(0, nrow = length(t), ncol = 9)
  lambda = matrix(0, nrow = length(t), ncol = 9)
  v1 = guess_v1
  v2 = guess_v2
  # final time adjoints
  lambda_init = rep(0,8);
  names(lambda_init) = paste0("lambda",1:8)
  # implement optimization
  oc = oc_optim(v1, v2, x, lambda, 
                IC, lambda_init, 
                bounds, delta, ode_fn, adj_fn, 
                t, params, oc_params)
  return(oc)
}

# implement loop for optimization
oc_optim = function(v1, v2, x, lambda, # initial guesses
                    IC, lambda_init, # state ICs & final time adjoints
                    bounds, delta, ode_fn, adj_fn,
                    t, params, oc_params){
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
      v1_interp <- approxfun(t, v1, rule = 2)
      v2_interp <- approxfun(t, v2, rule = 2)
      # solve states
      x <- ode(y = IC, times = t, func = ode_fn, parms = params, 
               v1_interp = v1_interp, v2_interp = v2_interp)
      # define interpolating functions for x (states)
      x_interp <- lapply(2:ncol(x), function(i){approxfun(x[,c(1,i)], rule = 2)})
      # solve adjoint equations (backwards)
      lambda <- ode(y = lambda_init, times = rev(t), func = adj_fn, 
                    parms = params, oc_params = oc_params, 
                    v1_interp = v1_interp, v2_interp = v2_interp, x_interp = x_interp, x = x)
      lambda <- lambda[nrow(lambda):1,]
      # calculate v1* and v2*
      temp_v1 <- ((lambda[,"lambda1"] - oc_params["C1"] - lambda[,"lambda3"])*x[,"S1"])/
        (2*oc_params["epsilon1"])
      temp_v2 <- ((lambda[,"lambda5"] - oc_params["C2"] - lambda[,"lambda7"])*x[,"S2"])/
        (2*oc_params["epsilon2"])
      # include bounds
      v1 = pmin(M1, pmax(0, temp_v1))
      v2 = pmin(M2, pmax(0, temp_v2))
      # update control
      v1 = 0.5*(v1 + oldv1)
      v2 = 0.5*(v2 + oldv2)
      # recalculate test
      test <- min(delta*norm_oc(c(v1,v2))-norm_oc(c(oldv1,oldv2)-c(v1,v2)),
                  delta*norm_oc(x[,-1])-norm_oc(oldx[,-1]-x[,-1]),
                  delta*norm_oc(lambda[,-1])-norm_oc(oldlambda[,-1]-lambda[,-1]))
      print(counter)
      print(test)
      counter <- counter + 1
    }
    return(list(x = x, lambda = lambda, v1 = v1, v2 = v2))
  })
}

# define norm(X,1) command from matlab
norm_oc <- function(x){sum(abs(x))}

# define cost function (j values)
j_integrand <- function(params, optim_states){
  with(as.list(c(params, optim_states)),{
    dj = b1 * (beta_I1*S1*I1 + beta_W1*S1*W1) + C1*v1*S1 + epsilon1*(v1^2) +
      b2 * (beta_I2*S2*I2 + beta_W2*S2*W2) + C2*v2*S2 + epsilon2*(v2^2)
    return(dj)
  })
}

calc_j <- function(params, optim_states, integrand_fn, 
                   lower_lim, upper_lim, step_size){
  # use trapezoid rule to integrate j
  steps <- seq(lower_lim, upper_lim, step_size)
  area = 0
  for(i in 1:(length(steps)-1)){
    a = steps[i]
    b = steps[i+1]
    j_a = integrand_fn(params, optim_states[i,])
    j_b = integrand_fn(params, optim_states[i+1,])
    area = area + (b-a)*(1/2)*(j_a+j_b)
  }
  return(area)
}

