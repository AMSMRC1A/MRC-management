# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)

source("code/cholera/CholeraSIRW_ODE.R")
source("code/cholera/Cholera_params.R")

# time varying vaccination
v1 = data.frame(times = t, v1 = rep(0,length(t)))
v1_interp <- approxfun(v1, rule = 2)
v2 = data.frame(times = t, v2 = rep(0,length(t)))
v2_interp <- approxfun(v2, rule = 2)

# solve ODE
out <- ode(y = IC, times = t, func = chol, parms = params, v1_interp = v1_interp, v2_interp = v2_interp)

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