## R resource
# http://desolve.r-forge.r-project.org/user2014/tutorial.pdf (see "Forcing" section)
# https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf

# initial guesses, conditions
# while
#   solve states 
#     interpolate control
#   solve adjoint
#     interpolate states and control over time

# initial guesses - controls
v1 = data.frame(times = times, v1 = rep(0,length(times)))
v1_interp <- approxfun(v1, rule = 2)
v1 = v1[,2]
v2 = data.frame(times = times, v2 = rep(0,length(times)))
v2_interp <- approxfun(v2, rule = 2)
v2 = v2[,2]
# adjoints
x = IC
lambda = rep(0,8)
names(lambda) = paste0("lambda",1:8)
# bounds
M1 = 0.015

# setup OC
delta = 0.001
test = -1
oc_params <- c(b1 = 1, b2 = 1,
               C1 = 1, C2 = 1, 
               epsilon1  = 1, epsilon2 = 1)

# define norm(X,1) command from matlab
norm <- function(x){sum(abs(x))}

counter <- 1

while(test < 0 & counter < 50){
  counter <- counter + 1
  # set previous control, state, and adjoint 
  oldv1 = v1
  oldv2 = v2
  oldx = x
  oldlambda = lambda
  
  # interpolate v
  v1_interp <- approxfun(v1, rule = 2)
  v2_interp <- approxfun(v2, rule = 2)
  
  # solve states
  x <- ode(y = x, times = t, func = chol, parms = params)
  
  # interpolate x (state)
  x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
  
  # solve adjoint
  lambda <- ode(y = lambda, times = t, func = adj, 
                parms = params, oc_params = oc_params,
                method = "bdf") # ?? implements backward differentiation?
  
  # calculate v1* and v2*
  temp_v1 <- ((lambda[,"lambda1"]- oc_params["C1"] - lambda[,"lambda3"])*x[,"S1"])/
    (2*oc_params["epsilon1"])
  temp_v2 <- ((lambda[,"lambda5"]- oc_params["C2"] - lambda[,"lambda7"])*x[,"S2"])/
    (2*oc_params["epsilon2"])
  v1 <- min(M1, max(0, temp_v1))
  v1 <- 0.5*(v1 + oldv1)
  
  # recalculate test
  test <- min(delta*norm(v1)-norm(oldv1-v1),
              delta*norm(x)-norm(oldy-x),
              delta*norm(lambda)-norm(oldlambda-lambda))
  
}


