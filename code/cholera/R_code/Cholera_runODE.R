# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)

source("CholeraSIRW_ODE.R")
source("Cholera_params.R")

#### SINGLE PARAMETER SETS (NO OPTIMAL CONTROL) ####

# time varying vaccination
v1 = data.frame(times = times, v1 = rep(0,length(times)))
v1_interp <- approxfun(v1, rule = 2)
v2 = data.frame(times = times, v2 = rep(0,length(times)))
v2_interp <- approxfun(v2, rule = 2)

# solve ODE
out <- ode(y = IC, times = times, func = chol, parms = params, v1_interp = v1_interp, v2_interp = v2_interp)

# reformat output for plotting
out <- as.data.frame(out)
out <- melt(out, id = c("time"))
out$compartment = substr(out$variable,1,1)
out$compartment = factor(out$compartment, levels = c("S", "I", "R", "W", "c", "V"))
out$patch = substr(out$variable, 2,2)

# plot output
ggplot(data = out, aes(x = time, y = value, color = patch))+
  geom_line(lwd=2)+
  facet_wrap(vars(compartment), scales = "free")+
  theme_half_open(12) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18))


#### MULTIPLE PARAMETER SETS (NO OPTIMAL CONTROL) ####

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(4)

# define parameters to test
test_params <- expand.grid(mu1 = c(0,0.01), # note changed this
                           beta_I1 = 0,#c(0,0.0000264), 
                           beta_W1 = 0.000121, #*c(1/10,1,10), going to assume this for now based off of Kyle's fig
                           m1 = 0.05*c(0,.5,1),
                           n1 = c(0,0.0001), 
                           gamma1 = 0.25,
                           delta1 = 5E-4, 
                           xi1 = 7.56E-3,
                           nu1 = 7.56E-3,
                           rho1 = c(0,0.025),
                           v1 = 0,
                           #mu2 = c(0,0.00001), 
                           beta_I2 = 0,# c(0,0.0000264), 
                           beta_W2 = 0.000121, #*c(1/10,1,10), going to assume this for now based off of Kyle's fig
                           m2 = 0.05*c(0,.5,1),
                           n2 = c(0,0.0001), 
                           gamma2 = 0.25,
                           delta2 = 5E-4, 
                           xi2 = 7.56E-3,
                           nu2 = 7.56E-3,
                           rho2 = c(0,0.025),
                           v2 = 0)
test_params$mu2 = test_params$mu1 # assume birth rates are equal right now
test_params$sim = 1:nrow(test_params)


new_times <- seq(0,100,0.1) # shorten time vector for memory

start <- Sys.time()
out <- foreach (i=1:nrow(test_params)) %dopar% { # can update so we don't have to include all params in test_params df
  ode(y = IC, times = new_times, func = chol, parms = test_params[i,])
}
end <- Sys.time()
end-start


# collect I1 and I2 for all cases
l <- length(new_times)
plot_df <- lapply(1:length(out), function(i){cbind(sim = rep(i, l), out[[i]][,c("time", "I1", "I2")])})
plot_df <- do.call(rbind, plot_df)
# add parameters of interest
plot_df <- left_join(test_params %>% select(mu1, m1, n1, rho1, m2, n2, rho2, sim),
                     as.data.frame(plot_df))
# transform to long for plotting
plot_df <- melt(plot_df, c("mu1", "m1", "n1", "rho1", "m2", "n2", "rho2", "sim", "time"), value.name = "count", variable.name = "class")

plts <- list()
for(i in unique(plot_df$rho2)){
  for(j in unique(plot_df$rho2))
  plts[[paste0(i,j)]] <- ggplot(data = plot_df %>% filter(rho1 == i, rho2 == j), aes(x = time, linetype = as.factor(mu1)))+
    geom_line(aes(y = count, color = class)) +
    ggtitle(paste0("rho1 =", i, ", rho2 = ", j), subtitle = "Assumptions: beta_I1 = beta_I2 = 0; beta_W1 = beta_W2 = 0.000121; mu1 = mu2")+
    facet_grid(cols = vars(m1, n1), rows = vars(m2, n2), labeller = "label_both")+
    scale_y_continuous(limits = range(plot_df$count)) +
    theme_classic() +
    theme(legend.position = "bottom")
}

pdf("figures/epicurves/epicurves_varying-mu-rho-m-n.pdf")
invisible(lapply(plts, print))
dev.off()




