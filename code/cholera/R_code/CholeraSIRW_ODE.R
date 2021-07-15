# model equations
chol <- function(times, y, params, v1_interp, v2_interp){
  with(as.list(c(y, params)),{
    v1 = v1_interp(times)
    v2 = v2_interp(times)
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