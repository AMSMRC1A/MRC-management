# define parameters
params <- c(mu1 = 0, mu2 = 0,                      # natural birth/death rate 1E-4
            beta_I1 = 2.14E-5, beta_I2 = 2.14E-5,  # transmission rate from people
            beta_W1 = 1.01E-5, beta_W2 = 1.01E-5,  # transmission rate from water
            v1 = 0, v2 = 0,                        # vaccination rate
            m1 = 0.05, m2 = 0.05,                  # movement rate (non-infected)
            n1 = 0, n2 = 0,                        # movement rate (infected)
            gamma1 = 0.25, gamma2 = 0.25,          # recovery rate
            delta1 = 5E-4, delta2 = 5E-4,          # disease induced mortality
            xi1 = 7.56E-3, xi2 = 7.56E-3,          # pathogen survival rate in water
            nu1 = 7.56E-3, nu2 = 7.56E-3,          # pathogen clearance rate in water
            rho1 = 0.05, rho2 = 0.05)#,            # pathogen movement rate in water
#b1 = 1, b2 = 1,                       # control cost for new cases
#C1 = 10, C2 = 10,                     # control cost for vaccination
#epsilon1 = 1, epsilon2 = 1)           # control quadratic term 

# define time (days)
t <- seq(0,200,0.05)



# define ICs
IC <- c(S1 = 10000-100, S2 = 10000-10,
        I1 = 100,    I2 = 10, 
        R1 = 0,      R2 = 0,
        W1 = 100,    W2 = 10)#,
# c1 = 0,      c2 = 0, # added c to count new cases
#  V1 = 0,      V2 = 0) # added V to count vaccinations

