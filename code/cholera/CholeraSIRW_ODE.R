# librarires
library(deSolve)
library(reshape2)

# model equations
chol <- function(t, y, params){
  with(as.list(c(y, params)),{
    dS1 <- mu1*(S1 + I1 + R1) - beta_I1*S1*I1 - beta_W1*S1*W1 - (mu1 + v1)*S1 - m1*S1 + m2*S2
    dS2 <- mu2*(S2 + I2 + R2) - beta_I2*S2*I2 - beta_W2*S2*W2 - (mu2 + v2)*S2 + m1*S1 - m2*S2
    dI1 <- beta_I1*S1*I1 + beta_W1*S1*W1 - (gamma1 + mu1 + delta1)*I1 - n1*I1 + n2*I2
    dI2 <- beta_I2*S2*I2 + beta_W2*S2*W2 - (gamma2 + mu2 + delta2)*I2 + n1*I1 - n2*I2
    dR1 <- gamma1*I1 - mu1*R1 + v1*S1 - m1*R1 + m2*R2
    dR2 <- gamma2*I2 - mu2*R2 + v2*S2 + m1*R1 - m2*R2
    dW1 <- xi1*I1 - v1*W1 - rho1*W1 + rho2*W2
    dW2 <- xi2*I2 - v2*W2 + rho1*W1 - rho2*W2
    ret <- c(dS1, dS2, dI1, dI2, dR1, dR2, dW1, dW2)
    return(list(ret))
  })
}

# define parameters
params <- c(mu1 = 1E-4, mu2 = 1E-4,                # natural birth/death rate
            beta_I1 = 2.64E-5, beta_I2 = 2.64E-5,  # transmission rate from people
            beta_W1 = 1.21E-4, beta_W2 = 1.21E-4,  # transmission rate from water
            v1 = 0, v2 = 0,                        # vaccination rate
            m1 = 0.1, m2 = 0.1,                    # movement rate (non-infected)
            n1 = 0, n2 = 0,                        # movement rate (infected)
            gamma1 = 0.25, gamma2 = 0.25,          # recovery rate
            delta1 = 5E-4, delta2 = 5E-4,          # disease induced mortality
            xi1 = 7.56E-3, xi2 = 7.56E-3,          # pathogen survival rate in water
            rho1 = 0.1, rho2 = 0.2)                # pathogen movement rate in water

# define ICs
IC <- c(S1 = 5000, S2 = 5000,
        I1 = 0,    I2 = 10, 
        R1 = 0,    R2 = 0,
        W1 = 0,    W2 = 25)

# define time (days)
t <- seq(0,120,0.01)

# solve ODE
out <- ode(y = IC, times = t, func = chol, parms = params)
out <- as.data.frame(out)
out <- melt(out, id = c("time"))
out$compartment = substr(out$variable,1,1)
out$compartment = factor(out$compartment, levels = c("S", "I", "R", "W"))
out$patch = substr(out$variable, 2,2)

# plot output
ggplot(data = out, aes(x = time, y = value, color = patch))+
  geom_line()+
  facet_wrap(vars(compartment), scales = "free")+
  theme_bw()

