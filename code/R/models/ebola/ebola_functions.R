# ODE model equations-----------------------------------------------------------

## Ebola_OC2.tex 
 
#' Ebola ODE model
#' 
#' defines 2-patch SEIHR model for ebola, with vaccination (v1, v2).
#' To be passed to \code{ode()}
#'  
#' @param t current time step
#' @param y vector of current conditions 
#' @param params vector of model parameters
#' @param interp_controls list of functions to interpolate controls 
#' (for time-varying controls), if NA constant control rates are assumed; names 
#' in list should correspond to variable name
ode_ebola <- function(t, y, params, interp_controls = NA) {
  with(as.list(c(y, params)), {
    # define population size
    N1<-S1+E1+I1+H1+R1
    N2<-S2+E2+I2+H2+R2
    # if controls are time-varying, use interpolation functions to generate
    # control rates over time
    if(all(!is.na(interp_controls))){
      # loop over all time-varying controls
      for(i in names(interp_controls)){
        assign(i, interp_controls[[i]](t))
      }
    }
    ## ODE model
    # patch 1
    dS1 <- mu1*N1 - betaI1*S1*I1 - betaD1*S1*D1 - mu1*S1 - v1*S1 - m1*S1 + m2*S2
    dE1 <- betaI1*S1*I1 + betaD1*S1*D1 - mu1*E1 - alpha1*E1 - m1*E1 + m2*E2
    dI1 <- alpha1*E1 - (mu1 + gammaI1 + (1 + u1)*phi1 + deltaI1 + n1)*I1 + n2*I2
    dH1 <- (1 + u1)*phi1*I1 - (gammaH1 + deltaH1 + mu1)*H1
    dD1 <- deltaI1*I1 - xi1*D1
    dR1 <- gammaI1*I1 + gammaH1*H1 + v1*S1 - mu1*R1 - m1*R1 + m2*R2
    # patch 2
    dS2 <- mu2*N2 - betaI2*S2*I2 - betaD2*S2*D2 - mu2*S2 - v2*S2 + m1*S1 - m2*S2
    dE2 <- betaI2*S2*I2 + betaD2*S2*D2 - mu2*E2 - alpha2*E2 + m1*E1 - m2*E2
    dI2 <- alpha2*E2 - (mu2 + gammaI2 + (1+u2)*phi2 + deltaI2)*I2 + n1*I1 - n2*I2
    dH2 <- (1 + u2)*phi2*I2 - (gammaH2 + deltaH2 + mu2)*H2
    dD2 <- deltaI2*I2 - xi2*D2
    dR2 <- gammaI2*I2 + gammaH2*H2 + v2*S2 - mu2*R2 + m1*R1 - m2*R2
    ret <- c(dS1,dE1,dI1,dH1,dD1,dR1,dS2,dE2,dI2,dH2,dD2,dR2)
    return(list(ret))
  })
}


# adjoint equations-------------------------------------------------------------

#' Ebola adjoint equations
#' 
#'  defines adjoint equations for ebola optimal control problem.
#'  To be passed to \code{ode()}
#' 
#'  @inheritParams ode_ebola
#'  @param x_interp list of functions to interpolate state variables, x
adjoint_ebola <- function(t, y, params, 
                            interp_controls = NA, 
                            x_interp, x) {
  # calculate state at time t using interpolated function
  state<- sapply(1:length(x_interp),function(i){x_interp[[i]](t)})
  names(state)<-c("S1","E1","I1","H1","D1","R1","S2","E2","I2","H2","D2","R2")
  # calculate control at time t using interpolated function
  if(all(!is.na(interp_controls))){
    # loop over all time-varying controls
    for(i in names(interp_controls)){
      assign(i, interp_controls[[i]](t))
    }
  }
  # adjoint equations
  with(as.list(c(y, params, state)), {
    # define population size
    N1<-S1+E1+I1+H1+R1
    N2<-S2+E2+I2+H2+R2
    # define adjoints
    # patch 1
    dlambda1 <- -( b1*(betaI1*I1 + betaD1*D1) + Cv1*v1 + 
                     lambda1*(mu1 -betaI1*I1 - betaD1*D1 - mu1 - v1 - m1)  + 
                     lambda2*(betaI1*I1 + betaD1*D1) + 
                     lambda6*v1 + 
                     lambda7*m1 )
    dlambda2 <- -( Cv1*v1  + 
                     lambda1*mu1 + 
                     lambda2*(-mu1 - alpha1 - m1) + 
                     lambda3*alpha1 + 
                     lambda8*m1 )
    dlambda3 <- -( b1*(betaI1*S1) + 
                     lambda1*(mu1 - betaI1*S1) + 
                     lambda2*(betaI1*S1) + 
                     lambda3*(-(mu1 + gammaI1 + (1+u1)*phi1 + deltaI1 + n1))  + 
                     lambda4*(1+u1)*phi1 + 
                     lambda5*deltaI1 + 
                     lambda6*gammaI1 + 
                     lambda9*(n1) )
    dlambda4 <- -(lambda1*(mu1) +
                    lambda4*(-(mu1 + gammaH1 + deltaH1)) + 
                    lambda6*(gammaH1))
    dlambda5 <- -( b1*(betaD1*S1) - 
                     lambda1*betaD1*S1 + 
                     lambda2*betaD1*S1 - 
                     lambda5*xi1)
    dlambda6 <- -( lambda1*mu1 + 
                     lambda6*(-mu1 - m1) + 
                     lambda12*m1)
    # patch 2
    dlambda7 <- -( b2*(betaI2*I2 + betaD2*D2) + Cv2*v2 + 
                     lambda7*(mu2 -betaI2*I2 - betaD2*D2 - mu2 - v2 - m2)  + 
                     lambda8*(betaI2*I2 + betaD2*D2) + 
                     lambda12*v2 + 
                     lambda1*m2)
    dlambda8 <- -( Cv2*v2  + 
                     lambda7*mu2 + 
                     lambda8*(-mu2 - alpha2 - m2) + 
                     lambda9*alpha2 + 
                     lambda2*m2)
    dlambda9 <- - (b2*betaI2*S2 + 
                     lambda7*(mu2 - betaI2*S2) + 
                     lambda8*(betaI2*S2) + 
                     lambda9*(-(mu2 + gammaI2 + (1+u2)*phi2 + deltaI2 + n2)) + 
                     lambda10*(1+u2)*phi2 + 
                     lambda11*deltaI2 + 
                     lambda12*gammaI2 + 
                     lambda3*n2)
    dlambda10 <- -(lambda7*mu2 + 
                     lambda10*(-(mu2 + gammaH2 + deltaH2)) + 
                     lambda12*gammaH2)
    dlambda11 <- -( b2*betaD2*S2 - 
                      lambda7*betaD2*S2 + 
                      lambda8*betaD2*S2 - 
                      lambda11*xi2 )
    dlambda12 <- -( lambda7*mu2 + 
                      lambda12*(-mu2 - m2) + 
                      lambda6*m2)
    ret <- c(
      dlambda1, dlambda2, dlambda3, dlambda4, dlambda5, dlambda6,
      dlambda7, dlambda8, dlambda9, dlambda10, dlambda11, dlambda12
    )
    return(list(ret))
  })
}

# optimal control solutions ----------------------------------------------------
 
#' calculate optimal control characterization
#' 
#' sub-function used in \code{update_optimal_solution()}
#' 
#' @param params vector of model parameters
#' @param lambda matrix of optimal lambda values (over time)
#' @param x matrix of optimal states (over time)
#' @param control_type character to define the type of control being implemented;
#' either \code{"uniform"} for the same control being applied in both patches 
#' or \code{"unique"} where control can vary across patches 
#' 
#' @return list of optimal control characterizations (vectors)
optimal_controls_ebola <- function(params, lambda, x, control_type){
  with(as.list(c(params, lambda, x)), {
    params <- as_tibble(as.list(params))
    if(control_type == "uniform")
      temp_controls <- list(
        v1 = (-Cv1*(S1+E1)+lambda1*S1-lambda6*S1 -
               Cv2*(S2+E2)+lambda7*S2-lambda12*S2)/(2*(epsilonV1+epsilonV2)),
        v2 = (-Cv1*(S1+E1)+lambda1*S1-lambda6*S1 -
               Cv2*(S2+E2)+lambda7*S2-lambda12*S2)/(2*(epsilonV1+epsilonV2)),
        u1 = (-Cu1*phi1*I1 + lambda3*phi1*I1 - lambda4*phi1*I1 + 
            -Cu2*phi1*I2 + lambda9*phi2*I2 - lambda10*phi2*I2)/
          (2*(epsilonU1+epsilonU2)),
        u2 = (-Cu1*phi1*I1 + lambda3*phi1*I1 - lambda4*phi1*I1 + 
                -Cu2*phi1*I2 + lambda9*phi2*I2 - lambda10*phi2*I2)/
          (2*(epsilonU1+epsilonU2))
      )
    else if(control_type == "unique"){
      temp_controls <- list(
        v1 = (-Cv1*(S1+E1)+lambda1*S1-lambda6*S1)/(2*epsilonV1),
        v2 = (-Cv2*(S2+E2)+lambda7*S2-lambda12*S2)/(2*epsilonV2),
        u1 = (-Cu1*phi1*I1 + lambda3*phi1*I1 - lambda4*phi1*I1)/(2*epsilonU1),
        u2 = (-Cu2*phi2*I2 + lambda9*phi2*I2 - lambda10*phi2*I2)/(2*epsilonU2)
      )
    }
    return(temp_controls)
  })
}

# EH: is there a better way to do this?
calc_test_ebola <- function(tol, controls, x, lambda, old_vals) {
  with(as.list(c(controls, x, lambda, old_vals)), {
    test_v <- tol * norm_oc(c(v1, v2)) - norm_oc(c(oldv1, oldv2) - c(v1, v2))
    test_u <- tol * norm_oc(c(u1, u2)) - norm_oc(c(oldu1, oldu2) - c(u1, u2))
    test_x <- tol * norm_oc(x[, -1]) - norm_oc(oldx[, -1] - x[, -1])
    test_lambda <- tol * norm_oc(lambda[, -1]) - norm_oc(oldlambda[, -1] - lambda[, -1])
    return(min(test_v, test_u, test_x, test_lambda))
  })
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
calc_j_ebola <- function(times, optim_states, params) {
  x <- times
  j_ints <- list(
    case1 = expression(b1*(betaI1*S1 + betaD1*S1*D1)),
    case2 = expression(b2*(betaI2*S2*I2 + betaD2*S2*D2)),
    vacc1 = expression(Cv1*v1*(S1+E1) + epsilonV1*(v1^2)),
    vacc2 = expression(Cv2*v2*(S1+E1) + epsilonV2*(v2^2)),
    sani1 = expression(Cu1*u1*phi1*I1 + epsilonU1*(u1^2)),
    sani2 = expression(Cu2*u2*phi2*I2 + epsilonU2*(u2^2))
  )
  j_vals <- lapply(j_ints, function(x) {
    apply(optim_states, 1, eval_j_integrand, params = params, integrand = x)
  })
  return(data.frame(
    j_case1 = trapz(x, j_vals[[1]]),
    j_case2 = trapz(x, j_vals[[2]]),
    j_vacc1 = trapz(x, j_vals[[3]]),
    j_vacc2 = trapz(x, j_vals[[4]]),
    j_sani1 = trapz(x, j_vals[[3]]),
    j_sani2 = trapz(x, j_vals[[4]])
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

# utilities --------------------------------------------------------------------

#' define norm(X,1) command from matlab
#' 
#' @param x vector
norm_oc <- function(x) {
  sum(abs(x))
}