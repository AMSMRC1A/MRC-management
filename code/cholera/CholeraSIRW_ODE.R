# libraries
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)

# model equations
chol <- function(t, y, params){
  with(as.list(c(y, params)),{
    v1 = v1_interp(t)
    v2 = v2_interp(t)
    dS1 <- mu1*(S1 + I1 + R1) - beta_I1*S1*I1 - beta_W1*S1*W1 - (mu1 + v1)*S1 - m1*S1 + m2*S2
    dS2 <- mu2*(S2 + I2 + R2) - beta_I2*S2*I2 - beta_W2*S2*W2 - (mu2 + v2)*S2 + m1*S1 - m2*S2
    dI1 <- beta_I1*S1*I1 + beta_W1*S1*W1 - (gamma1 + mu1 + delta1)*I1 - n1*I1 + n2*I2
    dI2 <- beta_I2*S2*I2 + beta_W2*S2*W2 - (gamma2 + mu2 + delta2)*I2 + n1*I1 - n2*I2
    dR1 <- gamma1*I1 - mu1*R1 + v1*S1 - m1*R1 + m2*R2
    dR2 <- gamma2*I2 - mu2*R2 + v2*S2 + m1*R1 - m2*R2
    dW1 <- xi1*I1 - nu1*W1 - rho1*W1
    dW2 <- xi2*I2 - nu2*W2 + rho1*W1 - rho2*W2
    #dc1 <- beta_I1*S1*I1 +  beta_W1*S1*W1 # cumulative infections patch 1
    #dc2 <- beta_I2*S2*I2 + beta_W2*S2*W2  # cumulative infections patch 2
    #dV1 <- v1*S1                          # cumulative vaccinations patch 1      
    #dV2 <- v2*S2                          # cumulative vaccinations patch 2      
    ret <- c(dS1, dS2, dI1, dI2, dR1, dR2, dW1, dW2)#, dc1, dc2, dV1, dV2)
    return(list(ret))
  })
}

# define parameters
params <- c(mu1 = 0, mu2 = 0,                      # natural birth/death rate 1E-4
            beta_I1 = 2.14E-5, beta_I2 = 2.14E-5,  # transmission rate from people
            beta_W1 = 1.01E-5, beta_W2 = 1.01E-5,  # transmission rate from water
            v1 = 0, v2 = 0,                        # vaccination rate
            m1 = 0.1, m2 = 0.1,                    # movement rate (non-infected)
            n1 = 0, n2 = 0,                        # movement rate (infected)
            gamma1 = 0.25, gamma2 = 0.25,          # recovery rate
            delta1 = 5E-4, delta2 = 5E-4,          # disease induced mortality
            xi1 = 7.56E-3, xi2 = 7.56E-3,          # pathogen survival rate in water
            nu1 = 7.56E-3, nu2 = 7.56E-3,          # pathogen clearance rate in water
            rho1 = 0.05, rho2 = 0.05)#,              # pathogen movement rate in water
            #b1 = 1, b2 = 1,                        # control cost for new cases
            #C1 = 10, C2 = 10,                      # control cost for vaccination
            #epsilon1 = 1, epsilon2 = 1)            # control quadratic term 

# define time (days)
t <- seq(0,200,0.01)

# time varying vaccination
v1 = data.frame(times = t, v1 = rep(0,length(t)))
v1_interp <- approxfun(v1, rule = 2)
v2 = data.frame(times = t, v2 = rep(0,length(t)))
v2_interp <- approxfun(v2, rule = 2)

# define ICs
IC <- c(S1 = 10000-10, S2 = 10000-10,
        I1 = 100,    I2 = 10, 
        R1 = 0,      R2 = 0,
        W1 = 300,    W2 = 10)#,
       # c1 = 0,      c2 = 0, # added c to count new cases
       #  V1 = 0,      V2 = 0) # added V to count vaccinations

# solve ODE
out <- ode(y = IC, times = t, func = chol, parms = params)
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

