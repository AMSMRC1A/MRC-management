# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)

# load files
source("CholeraSIRW_ODE.R")
source("Cholera_params.R")
source("Cholera_OCfunc.R")
source("Cholera_adjoints.R")

# initial guesses - controls
guess_v1 = rep(0,length(times))
guess_v2 = rep(0, length(times))

# adjoints
# x = matrix(0, nrow = length(times), ncol = 9)
# lambda = matrix(0, nrow = length(times), ncol = 9)
# lambda_init = rep(0,8)
# names(lambda_init) = paste0("lambda",1:8)
# bounds



# setup optimal control parameters
delta = 0.01
oc_params <- c(b1 = 1, b2 = 100,
               C1 = 0.125, C2 = 2.5, 
               epsilon1  = 100, epsilon2 = 1000)

# run optimization
oc = run_oc(guess_v1, guess_v2, IC, bounds, chol, adj,
              times, c(params, oc_params), delta)


# collect trajectories and controls
control_trajectories <- as.data.frame(oc$x)
control_trajectories$v1 = oc$v1
control_trajectories$v2 = oc$v2
control_trajectories <- melt(control_trajectories, id = c("time"))
control_trajectories$compartment = substr(control_trajectories$variable,1,1)
control_trajectories$compartment = factor(control_trajectories$compartment, levels = c("S", "I", "R", "W", "v"))
control_trajectories$patch = substr(control_trajectories$variable, 2,2)

# find no control ODE
out_ode <- ode(y = IC, times = times, func = chol, parms = params, 
           v1_interp = approxfun(data.frame(times = times, v1 = rep(0,length(times))), rule = 2), 
           v2_interp = approxfun(data.frame(times = times, v2 = rep(0,length(times))), rule = 2))
# find max control ODE
max_control_ode <- ode(y = IC, times = times, func = chol, parms = max_params, 
                       v1_interp = approxfun(data.frame(times = times, v1 = rep(bounds[[1]],length(times))), rule = 2), 
                       v2_interp = approxfun(data.frame(times = times, v2 = rep(bounds[[2]],length(times))), rule = 2))

# reformat output for plotting
out_ode <- as.data.frame(out_ode)
out <- melt(out_ode, id = c("time"))
out$compartment = substr(out$variable,1,1)
out$compartment = factor(out$compartment, levels = c("S", "I", "R", "W", "c", "V"))
out$patch = substr(out$variable, 2,2)

# plot trajectories and controls
control_plot <- control_trajectories %>% ggplot(aes(x = time, y = value, color = patch)) +
  geom_line(lwd=1.5)+
  #geom_line(data=out,aes(x = time, y = value, color = patch),linetype="dashed") +
  facet_wrap(vars(compartment), scales = "free") +
  theme_half_open(12) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18))
control_plot


# find j values
j <- calc_j(params = c(params, oc_params), 
            optim_states = cbind(oc$x, v1 = oc$v1, v2 = oc$v2), 
            integrand_fn = j_integrand, 
            lower_lim = min(times), upper_lim = max(times), step_size = 0.01)
j_no_control <- calc_j(params = c(params, oc_params), 
                       optim_states = cbind(out_ode, v1 = rep(0, nrow(out_ode)), v2 = rep(0, nrow(out_ode))), 
                       integrand_fn = j_integrand, 
                       lower_lim = min(times), upper_lim = max(times), step_size = 0.01)

j_max_control <- calc_j(params = c(max_params, oc_params), 
                       optim_states = cbind(max_control_ode, v1 = rep(bounds[[1]], nrow(max_control_ode)),
                                            v2 = rep(bounds[[2]], nrow(max_control_ode))), 
                       integrand_fn = j_integrand, 
                       lower_lim = min(times), upper_lim = max(times), step_size = 0.01)

print(paste("No control:", round(j_no_control,1),
            "Control: ", round(j,1),
            "Max control:", round(j_max_control,1),
            "rel change", round(j/j_no_control,3)))





#### test multiple oc parameter values ####

## create data.frameof parameters to test
# assume costs are equal in both patches
equal <- expand.grid(b1 = c(1,10), 
                     C1 = c(0.125, 0.625, 1.25), 
                     epsilon1 = c(1000, 100000)) # consider changing to 3E5
equal$b2 <- equal$b1
equal$C2 <- equal$C1
equal$epsilon2 <- equal$epsilon1
equal$scenario <- 1
# assume patch 2 has 5X cost of patch 1
p2greater <- expand.grid(b1 = c(1,10), 
                         C1 = c(0.125, 0.625, 1.25), 
                         epsilon1 = c(1000, 100000)) # consider changing to 3E5
p2greater$b2 <- 5*p2greater$b1
p2greater$C2 <- 5*p2greater$C1
p2greater$epsilon2 <- 5*p2greater$epsilon1
p2greater$scenario <- 2
# assume patch 1 has 5X cost of patch 2
p1greater <- expand.grid(b1 = c(1,10)*5, 
                         C1 = c(0.125, 0.625, 1.25)*5, 
                         epsilon1 = c(1000, 100000)*5) # consider changing to 3E5
p1greater$b2 <- 1/5*p1greater$b1
p1greater$C2 <- 1/5*p1greater$C1
p1greater$epsilon2 <- 1/5*p1greater$epsilon1
p1greater$scenario <- 3
# combine into single data.frame to run all combinations
test_params <- bind_rows(equal, p1greater, p2greater)
test_params$counter <- 1:nrow(test_params)

# calculate OC
start_time <- Sys.time()
mult_oc_params <- apply(test_params, 1, apply_oc, 
                        guess_v1 = guess_v1, guess_v2 = guess_v2, 
                        init_x = IC, bounds = bounds,
                        ode_fn = chol, adj_fn = adj,
                        times = times, params = c(params, oc_params), delta = delta)
mult_oc_params <- do.call(rbind, mult_oc_params)
end_time <- Sys.time()
end_time - start_time

# change to long for plotting
mult_oc_params <- melt(mult_oc_params %>% select(-counter), 
                       id = c("b1", "b2", "C1", "C2", "epsilon1", "epsilon2","j","scenario","time"))
# rename scenarios to be meaningful
mult_oc_params$scenario <- sapply(mult_oc_params$scenario, function(i){
  switch(i, "Equal cost", "Greater cost in P2", "Greater cost in P1")})

# add j values to plot
j_vals <- mult_oc_params %>% 
  group_by(b1, b2, C1, C2, epsilon1, epsilon2, scenario) %>% 
  summarise(j = unique(j)) %>% 
  mutate(label = paste0("j = ", round(j)), 
         x = 200, 
         y = max(bounds)*9/10)

# plot outcomes and save in a single file
plts <- list()
for(i in unique(mult_oc_params$scenario)){
  p <- ggplot(data = mult_oc_params %>% filter(scenario == i)) + 
    geom_line(aes(x = time, y = value, color = variable)) + 
    geom_text(data = j_vals %>% filter(scenario == i), aes(x = x, y = y, label = label), hjust = 1, vjust = 0) +
    ggtitle(i) +
    facet_grid(cols = vars(b1,b2, epsilon1, epsilon2), rows = vars(C1, C2), labeller = "label_both") +
    theme_half_open(12) +
    theme(legend.position = "bottom")
  plts[[i]] <- p
}

pdf("vary_cost_parameters.pdf")
invisible(lapply(plts, print))
dev.off()