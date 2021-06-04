require(deSolve)
require(tidyverse)
require(reshape2)
source("ebola_functions_mrc.R")
source("ebola_params_mrc.R")

times<-seq(0,200,by=.1)

test.sol<-ode(y=Y,
              times=times,
              func=ebola.ode,
              parms=parm,
              method = "ode45")
#plot(test.sol)


# initial guesses, conditions
# while
#   solve states 
#     interpolate control
#   solve adjoint
#     interpolate states and control over time

# initial guesses - controls
v1 <- data.frame(times = times, v1 = rep(0,length(times)))
v1_interp <- approxfun(v1, rule = 2)
v1 = v1[,2]
v2 = data.frame(times = times, v2 = rep(0,length(times)))
v2_interp <- approxfun(v2, rule = 2)
v2 = v2[,2]
# adjoints
x = matrix(0, nrow = length(times), ncol = 13)
lambda = matrix(0, nrow = length(times), ncol = 13)
lambda_init = rep(0,12)
names(lambda_init) = paste0("lambda",1:12)
# bounds
M1 = 0.1
M2 = 0.1

# setup OC

# define norm(X,1) command from matlab
norm <- function(x){sum(abs(x))}

delta<-.01
counter <- 0
test <- -1
while(test < 0 & counter < 10){
  counter <- counter + 1
  # set previous control, state, and adjoint 
  oldv1 <- v1
  oldv2 <- v2
  oldx <- x
  oldlambda <- lambda
  
  # interpolate v
  v1_interp <- approxfun(times, v1, rule = 2)
  v2_interp <- approxfun(times, v2, rule = 2)
  
  # solve states
  solx <- ode(y = Y, times = times, func = ebola.ode, parms = parm)
  
  # interpolate x (state)
  x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
  
  # solve adjoint
  lambda <- ode(y = lambda_init, times = rev(times), func = adjoints.ode, 
                parms = parm)
  
  lambda <- lambda[nrow(lambda):1,]
  
  # calculate v1* and v2*
  temp_v1 <- with(parm,(-C1*(solx[,"S1"]+solx[,"I1"])+lambda[,"lambda1"]*solx[,"S1"])/(2*epsilon1))
  temp_v2 <- with(parm,(-C1*(solx[,"S2"]+solx[,"I2"])+lambda[,"lambda7"]*solx[,"S2"])/(2*epsilon2))
  
  v1 <- pmin(M1, pmax(0, temp_v1))
  v1 <- 0.5*(v1 + oldv1)
  v2 <- pmin(M2, pmax(0, temp_v2))
  v2 <- 0.5*(v2 + oldv2)
  
  # recalculate test
  if(length(x) != length(oldx)){browser()}
  test <- min(delta*norm(c(v1,v2))-norm(c(oldv1,oldv2)-c(v1,v2)),
                        delta*norm(x[,-1])-norm(oldx[,-1]-x[,-1]),
                        delta*norm(lambda[,-1])-norm(oldlambda[,-1]-lambda[,-1]))
  print(counter)
  print(test)
}

J1=with(c(parm,data.frame(solx)),sum((b1*(betaI1*S1 + betaD1*S1*D1) + b2*(betaI2*S2*I2 + betaD2*S2*D2))*dt))
J2=with(c(parm,data.frame(solx)),sum(C1*v1*(S1)+C2*v2*S2+epsilon1*v1^2+epsilon2*v2^2))

plot(solx)

plot(solx[,"time"],v1)
plot(solx[,"time"],v2)
