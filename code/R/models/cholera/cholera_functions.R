
# ODE model equations-----------------------------------------------------------

#' cholera ODE model
#' 
#' defines 2-patch SIRW model for cholera, with two controls: 
#' vaccination (v1, v2) and sanitation (u1, u2). To be passed to \code{ode()}
#'  
#' @param t current time step
#' @param y vector of current conditions 
#' @param params vector of model parameters
#' @param interp_controls list of functions to interpolate controls 
#' (for time-varying controls), if NA constant control rates are assumed; names 
#' in list should correspond to variable name
ode_cholera <- function(t, y, params, interp_controls = NA) {
  with(as.list(c(y, params)), {
    # if controls are time-varying, use interpolation functions to generate
    # control rates over time
    if(any(is.function(interp_controls))){
      # loop over all time-varying controls
      for(i in names(interp_controls)){
        assign(i, interp_controls[[i]](t))
      }
    }
    # ODE model
    dS1 <- mu1*(S1+I1+R1) - beta_I1*S1*I1 - (1-u1)*beta_W1*S1*W1 - 
      (mu1 + v1)*S1 - m1*S1 + m2*S2
    dS2 <- mu2*(S2+I2+R2) - beta_I2*S2*I2 - (1-u2)*beta_W2*S2*W2 - 
      (mu2+v2)*S2 + m1*S1 - m2*S2
    dI1 <- beta_I1*S1*I1 + (1-u1)*beta_W1*S1*W1 - 
      (gamma1+mu1+delta1)*I1 - n1*I1 + n2*I2
    dI2 <- beta_I2*S2*I2 + (1-u2)*beta_W2*S2*W2 - 
      (gamma2+mu2+delta2)*I2 + n1*I1 - n2*I2
    dR1 <- gamma1*I1 - mu1*R1 + v1*S1 - m1*R1 + m2*R2
    dR2 <- gamma2*I2 - mu2*R2 + v2*S2 + m1*R1 - m2*R2
    dW1 <- xi1*I1 - nu1*W1 - rho1*W1
    dW2 <- xi2*I2 - nu2*W2 + rho1*W1 - rho2*W2
    ret <- c(dS1, dS2, dI1, dI2, dR1, dR2, dW1, dW2)
    return(list(ret))
  })
}


# adjoint equations-------------------------------------------------------------

#' Cholera adjoint equations
#' 
#'  defines adjoint equations for cholera optimal control problem.
#'  To be passed to \code{ode()}
#' 
#'  @inheritParams ode_cholera
#'  @param x_interp list of functions to interpolate state variables, x
adjoint_cholera <- function(t, y, params, 
                            interp_controls = NA, 
                            x_interp, x) {
  # calculate state at time t using interpolated function
  state <- c(
    S1 = x_interp[[1]](t), S2 = x_interp[[2]](t),
    I1 = x_interp[[3]](t), I2 = x_interp[[4]](t),
    R1 = x_interp[[5]](t), R2 = x_interp[[6]](t),
    W1 = x_interp[[7]](t), W2 = x_interp[[8]](t)
  )
  if(any(is.function(interp_controls))){
    # loop over all time-varying controls
    for(i in names(interp_controls)){
      assign(i, interp_controls[[i]](t))
    }
  }
  # adjoint equations
  with(as.list(c(y, params, state)), {
    dlambda1 <- -(b1*(beta_I1*I1 + (1-u1)*beta_W1*W1) + C1*v1 +
                    lambda1*(mu1 - beta_I1*I1 - (1-u1)*beta_W1*W1 - mu1 - v1 - m1) +
                    lambda2*(beta_I1*I1 + (1-u1)*beta_W1*W1) +
                    lambda3*v1 +
                    lambda5*m1)
    dlambda2 <- -(b1*beta_I1*S1 +
                    lambda1*(mu1 - beta_I1*S1) +
                    lambda2*(beta_I1*S1 - (gamma1 + mu1 + delta1 + n1)) +
                    lambda3*gamma1 +
                    lambda4*xi1 +
                    lambda6*n1)
    dlambda3 <- -(lambda1*mu1 -
                    lambda3*(mu1 + m1) +
                    lambda7*m1)
    dlambda4 <- -(b1*(1-u1)*beta_W1*S1 -
                    lambda1*(1-u1)*beta_W1*S1 +
                    lambda2*(1-u1)*beta_W1*S1 -
                    lambda4*(xi1 + rho1) +
                    lambda8*rho1)
    dlambda5 <- -(b2*(beta_I2*I2 + beta_W2*W2) + C2*v2 +
                    lambda1*m2 +
                    lambda5*(mu2 - beta_I2*I2 - (1-u2)*beta_W2*W2 - mu2 - v2 - m2) +
                    lambda6*(beta_I2*I2 + (1-u2)*beta_W2*W2) +
                    lambda7*v2)
    dlambda6 <- -(b2*beta_I2*S2 +
                    lambda2*n2 +
                    lambda5*(mu2 - beta_I2*S2) +
                    lambda6*(beta_I2*S2 - (gamma2 + mu2 + delta2 + n2)) +
                    lambda7*gamma2 +
                    lambda8*xi2)
    dlambda7 <- -(lambda3*m2 +
                    lambda5*mu2 -
                    lambda7*(mu2 + m2))
    dlambda8 <- -(b2*(1-u2)*beta_W2*S2 -
                    lambda5*(1-u2)*beta_W2*S2 +
                    lambda6*(1-u2)*beta_W2*S2 -
                    lambda8*(xi2 + rho2))
    ret <- c(
      dlambda1, dlambda2, dlambda3, dlambda4,
      dlambda5, dlambda6, dlambda7, dlambda8
    )
    return(list(ret))
  })
}

# optimal control solutions ----------------------------------------------------
 
#' calculate optimal vaccination strategy
#' 
#' sub-function used in \code{oc_optim()}
#' 
#' @param params vector of model parameters
#' @param lambda matrix of optimal lambda values (over time)
#' @param x matrix of optimal states (over time)
#' @param control_type character to define the type of control being implemented;
#' either \code{"uniform"} for the same control being applied in both patches 
#' or \code{"unique"} where control can vary across patches 
#' 
#' @return list of optimal vaccination in both patches (vectors)
opt_v_cholera <- function(params, lambda, x, control_type) {
  with(as.list(c(params, lambda, x)), {
    params <- as_tibble(as.list(params))
    # KD: !!! for consistency, change "lambda" to use same entry calls as params
    if (control_type == "uniform") {
      temp_v1 <- (((lambda1-C1-lambda3)*S1) + ((lambda5-C2-lambda7)*S2)) /
        (2*epsilon1 + 2*epsilon2)
      temp_v2 <- temp_v1
    }
    if (control_type == "unique") {
      temp_v1 <- ((lambda1-C1-lambda3)*S1) / (2*epsilon1)
      temp_v2 <- ((lambda5-C2-lambda7)*S2) / (2*epsilon2)
    }
    return(list(temp_v1 = temp_v1, temp_v2 = temp_v2))
  })
}

#' calculate optimal sanitation strategy
#' 
#' sub-function used in \code{oc_optim()}
#' 
#' @inheritParams opt_v_cholera
#' 
#' @return list of optimal sanitation in both patches (vectors)
opt_u_cholera <- function(params, lambda, x, control_type){
  with(as.list(c(params, lambda, x)), {
    params <- as_tibble(as.list(params))
    if (control_type == "uniform"){
      temp_u1 <- ((lambda2 + b1 - lambda1)*beta_W1*S1*W1 - D1 + 
                    (lambda6 + b2 - lambda5)*beta_W2*S2*W2 - D2) /
        (2*(eta1 + eta2))
      temp_u2 <- temp_u1
    }
    if (control_type == "unique"){
      temp_u1 <- ((beta_W1*S1*W1)*(b1 - lambda1 + lambda2) - D1) / (2*eta1)
      temp_u2 <- ((beta_W2*S2*W2)*(b2 - lambda5 + lambda6) - D2) / (2*eta2)
    }
    return(list(temp_u1 = temp_u1, temp_u2 = temp_u2))
  }
  )
}

# total cost -------------------------------------------------------------------

#' calculate total cost of a given strategy
#' 
#' sub-function used in \code{oc_optim()}
#' 
#' @param times vector of times over which optimal solution is defined
#' @param optim_states matrix of optimal states
#' @param params vector of model parameters
#' 
#' @return data.frame of each component of the total cost (cases, vaccination,
#' sanitation in patches 1 and 2)
calc_j_cholera <- function(times, optim_states, params) {
  x <- times
  j_ints <- list(
    case1 = expression(b1*(beta_I1*S1*I1 + (1-u1)*beta_W1*S1*W1)),
    case2 = expression(b2*(beta_I2*S2*I2 + (1-u2)*beta_W2*S2*W2)),
    vacc1 = expression(C1*v1*S1 + epsilon1*(v1^2)),
    vacc2 = expression(C2*v2*S2 + epsilon2*(v2^2)),
    sani1 = expression(D1*u1 + eta1*u1^2), 
    sani2 = expression(D2*u2 + eta2*u2^2)
  )
  j_vals <- lapply(j_ints, function(x) {
    apply(optim_states, 1, eval_j_integrand, params = params, integrand = x)
  })
  return(data.frame(
    j_case1 = trapz(x, j_vals[[1]]),
    j_case2 = trapz(x, j_vals[[2]]),
    j_vacc1 = trapz(x, j_vals[[3]]),
    j_vacc2 = trapz(x, j_vals[[4]]),
    j_sani1 = trapz(x, j_vals[[5]]),
    j_sani2 = trapz(x, j_vals[[6]])
  ))
}

#' helper function to calculate costs
#' 
#' EH NOTE: THIS IS A GENERAL HELPER FUNCTION, DOES IT BELONG HERE?
#' 
#' @inheritParams calc_j
#' @param integrand list of integrands to be evaluated
#' 
#' @return list of evaluated integrands
eval_j_integrand <- function(params, optim_states, integrand) {
  with(as.list(c(optim_states, params)), {
    eval(parse(text = integrand))
  })
}


# R0 calculator ----------------------------------------------------------------

#' calculate R0 for the cholera model
#' 
#' @param params vector of model parameters
#' @params N0_1 integer (or double?) of initial population size in patch 1
#' @params N0_2 integer (or double?) of initial population size in patch 2
#' 
#' @return double of R0 value
R0_cholera <- function(params, N0_1, N0_2) {
  with(as.list(params), {
    # calculate intermediate parameters
    lambda1 <- 1 / (gamma1 + mu1 + delta1 + n1)
    lambda2 <- 1 / (gamma2 + mu2 + delta2 + n2)
    sigma1 <- n1 * lambda1
    sigma2 <- n2 * lambda2
    tau <- 1 / (1 - (sigma1 * sigma2))
    eta1 <- 1 / (nu1 + rho1)
    eta2 <- 1 / (nu2 + rho2)
    omega1 <- rho1 / (nu1 + rho1)
    omega2 <- rho2 / (nu2 + rho2)
    # calculate RIi and RWi
    RI1 <- beta_I1 * N0_1 * lambda1
    RI2 <- beta_I2 * N0_2 * lambda2
    RW1 <- beta_W1 * N0_1 * xi1 * lambda1 * eta1
    RW2 <- beta_W2 * N0_2 * xi2 * lambda2 * eta2
    # calculate components of R0
    R11 <- tau * (RI1 + RW1)
    R12 <- sigma2 * tau * (RI1 + RW1)
    R21 <- sigma1 * tau * (RI2 + RW2) + tau * omega1 * ((xi1 * lambda1) / (xi2 * lambda2)) * RW2
    R22 <- tau * (RI2 + RW2) + tau * sigma2 * omega1 * ((xi1 * lambda1) / (xi2 * lambda2)) * RW2
    R0 <- (R11 + R22 + sqrt((R11 - R22)^2 + 4 * R12 * R21)) / 2
    return(R0)
  })
}

