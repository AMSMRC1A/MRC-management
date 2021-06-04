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

#source("/src/cholera/Cholera_adjoints.R")

# initial guesses - controls
v1 = data.frame(times = t, v1 = rep(0,length(t)))
v1_interp <- approxfun(v1, rule = 2)
v1 = v1[,2]
v2 = data.frame(times = t, v2 = rep(0,length(t)))
v2_interp <- approxfun(v2, rule = 2)
v2 = v2[,2]
# adjoints
x = matrix(0, nrow = length(t), ncol = 9)
lambda = matrix(0, nrow = length(t), ncol = 9)
lambda_init = rep(0,8)
names(lambda_init) = paste0("lambda",1:8)
# bounds
M1 = 0.015
M2 = 0.015

# setup optimal control parameters
delta = 0.01
test = -1
oc_params <- c(b1 = 1, b2 = 1,
               C1 = 0.125, C2 = 0.125, 
               epsilon1  = 10, epsilon2 = 10)

# define norm(X,1) command from matlab (!! over-writes R "norm" function)
norm <- function(x){sum(abs(x))}

counter <- 1
test = -1
while(test < 0 & counter < 50){
  counter <- counter + 1
  # set previous control, state, and adjoint 
  oldv1 = v1
  oldv2 = v2
  oldx = x
  oldlambda = lambda
  
  # define interpolating functions for v
  v1_interp <- approxfun(t, v1, rule = 2)
  v2_interp <- approxfun(t, v2, rule = 2)
  
  # solve states
  x <- ode(y = IC, times = t, func = chol, parms = params)
  
  # define interpolating functions for x (states)
  x_interp <- lapply(2:ncol(x), function(i){approxfun(x[,c(1,i)], rule = 2)})
  
  # solve adjoint equations (backwards)
  lambda <- ode(y = lambda_init, times = rev(t), func = adj, 
                parms = params, oc_params = oc_params)
  lambda <- lambda[nrow(lambda):1,]
  
  # calculate v1* and v2*
  temp_v1 <- ((lambda[,"lambda1"] - oc_params["C1"] - lambda[,"lambda3"])*x[,"S1"])/
    (2*oc_params["epsilon1"])
  temp_v2 <- ((lambda[,"lambda5"] - oc_params["C2"] - lambda[,"lambda7"])*x[,"S2"])/
    (2*oc_params["epsilon2"])
  
  v1 <- pmin(M1, pmax(0, temp_v1))
  v1 <- 0.5*(v1 + oldv1)
  v2 <- pmin(M2, pmax(0, temp_v2))
  v2 <- 0.5*(v2 + oldv2)
  
  # recalculate test
  test <- min(delta*norm(c(v1,v2))-norm(c(oldv1,oldv2)-c(v1,v2)),
              delta*norm(x[,-1])-norm(oldx[,-1]-x[,-1]),
              delta*norm(lambda[,-1])-norm(oldlambda[,-1]-lambda[,-1]))
  print(counter)
  print(test)
}

# Collect trajectories and controls
control_trajectories <- as.data.frame(x)
control_trajectories$v1 = v1
control_trajectories$v2 = v2
control_trajectories <- melt(control_trajectories, id = c("time"))
control_trajectories$compartment = substr(control_trajectories$variable,1,1)
control_trajectories$compartment = factor(control_trajectories$compartment, levels = c("S", "I", "R", "W", "v"))
control_trajectories$patch = substr(control_trajectories$variable, 2,2)


# Plot trajectories and controls
control_plot <- control_trajectories %>% ggplot(aes(x = time, y = value, color = patch)) +
  geom_line(lwd=1.5)+
  geom_line(data=out,aes(x = time, y = value, color = patch),linetype="dashed") +
  facet_wrap(vars(compartment), scales = "free") +
  theme_half_open(12) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18))

control_plot