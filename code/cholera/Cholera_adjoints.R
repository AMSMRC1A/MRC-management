# adjoints
adj <- function(t, y, params, oc_params){
  # calculate state at time t using interpolated function
  state = sapply(1:length(x_interp), function(i){x_interp[[i]](t)})
  names(state) = colnames(x)[-1]
  v1 = v1_interp(t)
  v2 = v2_interp(t)
  # adjoint equations
  with(as.list(c(y, params, oc_params, state)),{
    dlambda1 = -(b1*(beta_I1*I1 + beta_W1*W1) + C1*v1 +
                   lambda1*(mu1 - beta_I1*I1 - beta_W1*W1 - mu1 - v1 - m1)+
                   lambda2*(beta_I1*I1 + beta_W1*W1)+
                   lambda3*v1 + 
                   lambda5*m1)
    dlambda2 = -(b1*beta_I1*S1 + 
                   lambda1*(mu1 - beta_I1*S1)+
                   lambda2*(beta_I1*S1 - (gamma1 + mu1 + delta1 + n1))+
                   lambda3*gamma1 +
                   lambda4*xi1 +
                   lambda6*n1)
    dlambda3 = -(lambda1*mu1 -
                   lambda3*(mu1 + m1)+
                   lambda7*m1)
    dlambda4 = -(b1*beta_W1*S1 - 
                   lambda1*beta_W1*S1  + 
                   lambda2*beta_W1*S1 -
                   lambda4*(xi1 + rho1)+
                   lambda8*rho1)
    dlambda5 = -(b2*(beta_I2*I2 + beta_W2*W2) + C2*v2 +
                   lambda1*m2 +
                   lambda5*(mu2 - beta_I2*I2 - beta_W2*W2 - mu2 - v2 - m2)+
                   lambda6*(beta_I2*I2 + beta_W2*W2)+
                   lambda7*v2)
    dlambda6 = -(b2*beta_I2*S2 + 
                   lambda2*n2 +
                   lambda5*(mu2 - beta_I2*S2)+
                   lambda6*(beta_I2*S2 - (gamma2 + mu2 + delta2 + n2))+
                   lambda7*gamma2 +
                   lambda8*xi2)
    dlambda7 = -(lambda3*m2 + 
                   lambda5*mu2 -
                   lambda7*(mu2 + m2))
    dlambda8 = -(b2*beta_W2*S2 - 
                   lambda5*beta_W2*S2 + 
                   lambda6*beta_W2*S2 -
                   lambda8*(xi2 + rho2))
    ret <- c(dlambda1, dlambda2, dlambda3, dlambda4,
             dlambda5, dlambda6, dlambda7, dlambda8)
    return(list(ret))
  })
}

