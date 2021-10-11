# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)
library(pracma)

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(4) # update this 4 if you want to use more cores

# load files
source("CholeraSIRW_ODE.R")
source("Cholera_params.R")
source("Cholera_OCfunc.R")
source("Cholera_adjoints.R")

# initial guesses - controls
guess_v1 <- rep(0, length(times))
guess_v2 <- rep(0, length(times))

# adjoints
# x = matrix(0, nrow = length(times), ncol = 9)
# lambda = matrix(0, nrow = length(times), ncol = 9)
# lambda_init = rep(0,8)
# names(lambda_init) = paste0("lambda",1:8)
# bounds



# setup optimal control parameters
delta <- 0.01
oc_params <- c(
  b1 = 1, b2 = 1,
  C1 = 0.125, C2 = 0.125,
  epsilon1 = 10000, epsilon2 = 10000
)

# run optimization
oc <- run_oc(
  guess_v1, guess_v2, IC, bounds, chol, adj,
  times, as.list(c(params, oc_params)), delta, "unique"
)

# collect trajectories and controls
control_trajectories <- as.data.frame(oc$x)
control_trajectories$v1 <- oc$v1
control_trajectories$v2 <- oc$v2
control_trajectories <- melt(control_trajectories, id = c("time"))
control_trajectories$compartment <- substr(control_trajectories$variable, 1, 1)
control_trajectories$compartment <- factor(control_trajectories$compartment, levels = c("S", "I", "R", "W", "v"))
control_trajectories$patch <- substr(control_trajectories$variable, 2, 2)

# find no control ODE
out_ode <- ode(
  y = IC, times = times, func = chol, parms = params,
  v1_interp = approxfun(data.frame(times = times, v1 = rep(0, length(times))), rule = 2),
  v2_interp = approxfun(data.frame(times = times, v2 = rep(0, length(times))), rule = 2)
)
# find max control ODE
max_control_ode <- ode(
  y = IC, times = times, func = chol, parms = max_params,
  v1_interp = approxfun(data.frame(times = times, v1 = rep(bounds[[1]], length(times))), rule = 2),
  v2_interp = approxfun(data.frame(times = times, v2 = rep(bounds[[2]], length(times))), rule = 2)
)

# reformat output for plotting
out_ode <- as.data.frame(out_ode)
out <- melt(out_ode, id = c("time"))
out$compartment <- substr(out$variable, 1, 1)
out$compartment <- factor(out$compartment, levels = c("S", "I", "R", "W", "c", "V"))
out$patch <- substr(out$variable, 2, 2)

# plot trajectories and controls
control_plot <- control_trajectories %>% ggplot(aes(x = time, y = value, color = patch)) +
  geom_line(lwd = 1.5) +
  # geom_line(data=out,aes(x = time, y = value, color = patch),linetype="dashed") +
  facet_wrap(vars(compartment), scales = "free") +
  theme_half_open(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18)
  )
control_plot


# find j values
j <- calc_j(
  params = c(params, oc_params), times = times,
  optim_states = cbind(oc$x, v1 = oc$v1, v2 = oc$v2), #< length of this (=4001) must match seq(lower_lim, upper_lim, step_size) = 20001
  integrand_fn = j_integrand
) # ,
# lower_lim = min(times), upper_lim = max(times), step_size = 0.01)
j_no_control <- calc_j(
  params = c(params, oc_params), times = times,
  optim_states = cbind(out_ode, v1 = rep(0, nrow(out_ode)), v2 = rep(0, nrow(out_ode))),
  integrand_fn = j_integrand
) # ,
# lower_lim = min(times), upper_lim = max(times), step_size = 0.01)
j_max_control <- calc_j(
  params = c(max_params, oc_params), times = times,
  optim_states = cbind(max_control_ode,
    v1 = rep(bounds[[1]], nrow(max_control_ode)),
    v2 = rep(bounds[[2]], nrow(max_control_ode))
  ),
  integrand_fn = j_integrand
) # ,
# lower_lim = min(times), upper_lim = max(times), step_size = 0.01)

print(paste(
  "No control:", round(j_no_control, 1),
  "Control: ", round(j, 1),
  "Max control:", round(j_max_control, 1),
  "rel change", round(j / j_no_control, 3)
))
