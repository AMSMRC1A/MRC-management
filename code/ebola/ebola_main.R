require(deSolve)
require(tidyverse)
require(reshape2)
source("ebola_functions_mrc.R")
source("ebola_params_mrc.R")
(R0<-with(c(parm,Y),sum(S1+E1+I1+R1+H1)*alpha1/(alpha1 + mu1)*1/(gammaI1 + phi1 + deltaI1 + mu1)*(betaI1 + betaD1*deltaI1/xi1)))
##------------------------------------------
times<-seq(0,730,by=2)

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
parm<-c(parm,v1_interp=v1_interp,v2_interp=v2_interp)
# adjoints
x = matrix(0, nrow = length(times), ncol = 13)
lambda = matrix(0, nrow = length(times), ncol = 13)
lambda_init = rep(0,12)
names(lambda_init) = paste0("lambda",1:12)
# bounds
M1 = 0.01
M2 = 0.01

# setup OC

#force same vaccination rates
v1=v2

# define norm(X,1) command from matlab
##------------------------------------
norm_oc <- function(x){colSums(abs(x))}

delta<-.01
counter <- 0
test <- -1
solx<-ode(y=Y,
          times=times,
          func=ebola.ode,
          parms=parm,
          method = "ode45")
x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
parm<-c(parm,x_interp=x_interp)
while(test < 0 & counter < 200){
  counter <- counter + 1
  # set previous control, state, and adjoint 
  oldv1 <- v1
  oldv2 <- v2
  oldx <- solx
  oldlambda <- lambda
  
  # interpolate v
  parm$v1_interp <- approxfun(times, v1, rule = 2)
  parm$v2_interp <- approxfun(times, v2, rule = 2)
  
  
  # solve states
  solx <- ode(y = Y, times = times, func = ebola.ode, parms = parm)
  
  # interpolate x (state)
  parm$x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
  
  # solve adjoint
  lambda <- ode(y = lambda_init, times = rev(times), func = adjoints.ode, 
                parms = parm)
  
  lambda <- lambda[nrow(lambda):1,]
  
  # calculate v1* and v2*
  temp_v1 <- with(parm,(-C1*(solx[,"S1"]+solx[,"I1"])+lambda[,"lambda1"]*solx[,"S1"]-lambda[,"lambda6"]*solx[,"S1"]-
                          C2*(solx[,"S2"]+solx[,"I2"])+lambda[,"lambda7"]*solx[,"S2"]-lambda[,"lambda12"]*solx[,"S2"])/(2*epsilon1))
#  temp_v2 <- with(parm,()/(2*epsilon2))
  
  v1 <- pmin(M1, pmax(0, temp_v1))
  v1 <- 0.5*(v1 + oldv1)
  v2 <- pmin(M2, pmax(0, temp_v1))
  v2 <- 0.5*(v1 + oldv1)
  
  # recalculate test
  if(length(x) != length(oldx)){browser()}
  
  test <- min(delta*norm_oc(data.frame(v1,v2))-norm_oc(data.frame(oldv1,oldv2)-data.frame(v1,v2)),
                        delta*norm_oc(solx[,-1])-norm_oc(oldx[,-1]-solx[,-1]),
                        delta*norm_oc(lambda[,-1])-norm_oc(oldlambda[,-1]-lambda[,-1]))
  print(c(delta*norm_oc(data.frame(v1,v2))-norm_oc(data.frame(oldv1,oldv2)-data.frame(v1,v2)),
          delta*norm_oc(solx[,-1])-norm_oc(oldx[,-1]-solx[,-1]),
          delta*norm_oc(lambda[,-1])-norm_oc(oldlambda[,-1]-lambda[,-1])))
  print(counter)
  print(test)
}
dt<-.1
(J1=with(c(parm,data.frame(solx)),sum((b1*(betaI1*S1 + betaD1*S1*D1) + b2*(betaI2*S2*I2 + betaD2*S2*D2))*dt)))
(J2=with(c(parm,data.frame(solx)),sum(C1*v1*(S1)+C2*v2*S2+epsilon1*v1^2+epsilon2*v2^2)*dt))

optimal_sol<-solx%>%
  data.frame(v1=v1,v2=v2)%>%
  pivot_longer(cols=names(.)[2:dim(.)[2]],names_to="State",values_to="value")%>%
  mutate(patch=substr(State,2,2),State=substr(State,1,1))

optimal_sol$State<-optimal_sol$State%>%
  factor(levels=c("S","E","I","R","H","D","v"),ordered = T)

ggplot(optimal_sol,aes(x=time,y=value,color=patch))+
  geom_line()+
  facet_wrap(~State,scales = "free")+
  theme_bw()
