ebola.ode<-function(t,Y,p){
  with(as.list(c(Y,p)),{
    N1<-S1+E1+I1+H1+R1
    N2<-S2+E2+I2+H2+R2
    
    v1<-v1_interp(t)
    v2<-v2_interp(t)
    
    dS1 <- mu1*N1 - betaI1*S1*I1 - betaD1*S1*D1 - mu1*S1 - v1*S1 - m1*S1 + m2*S2
    dE1 <- betaI1*S1*I1 + betaD1*S1*D1 - mu1*E1 - alpha1*E1 - m1*E1 + m2*E2
    dI1 <- alpha1*E1 - (mu1 + gammaI1 + phi1 + deltaI1)*I1 - n1*I1 + n2*I2
    dH1 <- phi1*I1 - (gammaH1 + deltaH1 + mu1)*H1
    dD1 <- deltaI1*I1 - xi1*D1
    dR1 <- gammaI1*I1 + gammaH1*H1 + v1*S1 - mu1*R1 - m1*R1 + m2*R2
    
    dS2 <- mu2*N2 - betaI2*S2*I2 - betaD2*S2*D2 - mu2*S2 - v2*S2 + m1*S1 - m2*S2
    dE2 <- betaI2*S2*I2 + betaD2*S2*D2 - mu2*E2 - alpha2*E2 + m1*E1 - m2*E2
    dI2 <- alpha2*E2 - (mu2 + gammaI2 + phi2 + deltaI2)*I2 + n1*I1 - n2*I2
    dH2 <- phi2*I2 - (gammaH2 + deltaH2 + mu2)*H2
    dD2 <- deltaI2*I2 - xi2*D2
    dR2 <- gammaI2*I2 + gammaH2*H2 + v2*S2 - mu2*R2 + m1*R1 - m2*R2
    return(list(c(dS1,dE1,dI1,dH1,dD1,dR1,dS2,dE2,dI2,dH2,dD2,dR2)))
  })
}

adjoints.ode<-function(t,Y,p,original_inits){
  
  state<- sapply(1:length(p$x_interp),function(i){p$x_interp[[i]](t)})
  names(state)<-c("S1","E1","I1","H1","D1","R1","S2","E2","I2","H2","D2","R2")
  with(as.list(c(Y,p,state)),{
    
    N1<-S1+E1+I1+H1+R1
    N2<-S2+E2+I2+H2+R2
    
    v1<-v1_interp(t)
    v2<-v2_interp(t)
    
    dlambda1 <- -( b1*(betaI1*I1 + betaD1*D1) + C1*v1 + lambda1*(mu1 -betaI1*I1 - betaD1*D1 - mu1 - v1 - m1)  + lambda2*(betaI1*I1 + betaD1*D1) + lambda6*v1 + lambda7*m1 );
    dlambda2 <- -( C1*v1  + lambda1*mu1 + lambda2*(-mu1 - alpha1 - m1) + lambda3*alpha1 + lambda8*m1 );
    dlambda3 <- -( b1*(betaI1*S1) + lambda1*(mu1 - betaI1*S1) + lambda2*(betaI1*S1) + lambda3*(-(mu1 + gammaI1 + phi1 + deltaI1 + n1))  + lambda4*phi1 + lambda5*deltaI1 + lambda6*gammaI1 + lambda9*(n1) );
    dlambda4 <- -(lambda1*(mu1) + lambda4*(-(mu1 + gammaH1 + deltaH1)) + lambda6*(gammaH1));
    dlambda5 <--( b1*(betaD1*S1) - lambda1*betaD1*S1 + lambda2*betaD1*S1 - lambda5*xi1);
    dlambda6 <- -( lambda1*mu1 + lambda6*(-mu1 - m1) + lambda12*m1);
    
    # patch 2
    dlambda7 <- -( b2*(betaI2*I2 + betaD2*D2) + C2*v2 + lambda7*(mu2 -betaI2*I2 - betaD2*D2 - mu2 - v2 - m2)  + lambda8*(betaI2*I2 + betaD2*D2) + lambda12*v2 + lambda1*m2);
    dlambda8 <- -( C2*v2  + lambda7*mu2 + lambda8*(-mu2 - alpha2 - m2) + lambda9*alpha2 + lambda2*m2);
    dlambda9 <- - (b2*betaI2*S2 + lambda7*(mu2 - betaI2*S2) + lambda8*(betaI2*S2) + lambda9*(-(mu2 + gammaI2 + phi2 + deltaI2 + n2)) + lambda10*phi2 + lambda11*deltaI2 + lambda12*gammaI2 + lambda3*n2); 
    dlambda10 <- -(lambda7*mu2 + lambda10*(-(mu2 + gammaH2 + deltaH2)) + lambda12*gammaH2);
    dlambda11 <- -( b2*betaD2*S2 - lambda7*betaD2*S2 + lambda8*betaD2*S2 - lambda11*xi2 );
    dlambda12 <- -( lambda7*mu2 + lambda12*(-mu2 - m2) + lambda6*m2);
    return(list(c(dlambda1,dlambda2,dlambda3,dlambda4,dlambda5,dlambda6,
                  dlambda7,dlambda8,dlambda9,dlambda10,dlambda11,dlambda12)))
  })
}


##This function replicates what was previously in the main code.  It takes in
##initial conditions (inits), parameters (params) as they are currently specified
##in ebola_params_mrc.R, a vector of the maximum vacination rates M=c(M1,M2), and
##the time period for the optimal control problem (defaults to 2 years if missing)
##It returns the results in a named list with the simulated result using the 
##optimal control,including vaccination rates, (SolX) and the final costs 
##(J11 (Disease in patch 1), J12 (Disease in patch 2), J21 (vaccination in patch 1)
##, and J22 (vaccination in patch 2)),

ebola.optim<-function(inits,params,M,times=seq(0,730,by=1)){
  # initial guesses - controls
  v1 <- data.frame(times = times, v1 = rep(0,length(times)))
  v1_interp <- approxfun(v1, rule = 2)
  v1 = v1[,2]
  v2 = data.frame(times = times, v2 = rep(0,length(times)))
  v2_interp <- approxfun(v2, rule = 2)
  v2 = v2[,2]
  params<-c(params,v1_interp=v1_interp,v2_interp=v2_interp)
  # adjoints
  x = matrix(0, nrow = length(times), ncol = 13)
  lambda = matrix(0, nrow = length(times), ncol = 13)
  lambda_init = rep(0,12)
  names(lambda_init) = paste0("lambda",1:12)
  
  # define norm(X,1) command from matlab
  ##------------------------------------
  norm_oc <- function(x){colSums(abs(x))}
  
  delta<-.01
  counter <- 0
  test <- -1
  solx<-ode(y=Y,
            times=times,
            func=ebola.ode,
            parms=params,
            method = "ode45")
  x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
  params<-c(params,x_interp=x_interp)
  while(test < 0 & counter < 200){
    counter <- counter + 1
    # set previous control, state, and adjoint 
    oldv1 <- v1
    oldv2 <- v2
    oldx <- solx
    oldlambda <- lambda
    
    # interpolate v
    params$v1_interp <- approxfun(times, v1, rule = 2)
    params$v2_interp <- approxfun(times, v2, rule = 2)
    
    
    # solve states
      solx <- ode(y = inits, times = times, func = ebola.ode, parms = params)
    
    
    # interpolate x (state)
    params$x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
    
    
    # solve adjoint
    lambda <- ode(y = lambda_init, times = rev(times), func = adjoints.ode, 
                  parms = params)
    
    lambda <- lambda[nrow(lambda):1,]
    
    # calculate v1* and v2*
    temp_v1 <- with(params,(-C1*(solx[,"S1"]+solx[,"I1"])+lambda[,"lambda1"]*solx[,"S1"]-lambda[,"lambda6"]*solx[,"S1"])/(2*epsilon1))
    temp_v2 <- with(params,(-C2*(solx[,"S2"]+solx[,"I2"])+lambda[,"lambda7"]*solx[,"S2"]-lambda[,"lambda12"]*solx[,"S2"])/(2*epsilon2))
    
    v1 <- pmin(M[1], pmax(0, temp_v1))
    v1 <- 0.5*(v1 + oldv1)
    v2 <- pmin(M[2], pmax(0, temp_v2))
    v2 <- 0.5*(v2 + oldv2)
    
    # recalculate test
    if(length(x) != length(oldx)){browser()}
    
    test <- min(delta*norm_oc(data.frame(v1,v2))-norm_oc(data.frame(oldv1,oldv2)-data.frame(v1,v2)),
                delta*norm_oc(solx[,-1])-norm_oc(oldx[,-1]-solx[,-1]),
                delta*norm_oc(lambda[,-1])-norm_oc(oldlambda[,-1]-lambda[,-1]))
    print(counter)
    print(test)
    
  }
  dt<-.1
  J11=with(c(params,data.frame(solx)),sum((b1*(betaI1*S1 + betaD1*S1*D1) )*dt))
  J12=with(c(params,data.frame(solx)),sum((b2*(betaI2*S2*I2 + betaD2*S2*D2))*dt))
  J21=with(c(params,data.frame(solx)),sum(C1*v1*(S1)+epsilon1*v1^2)*dt)
  J22=with(c(params,data.frame(solx)),sum(C2*v2*S2+epsilon2*v2^2)*dt)
  solx<-data.frame(solx,v1,v2)
  returned<-list(solx,J11,J12,J21,J22)
  names(returned)<-c("SolX","J11","J12","J21","J22")
  return(returned)
}


##Optimal control for simulation where infection is introduced (I=10) in patch 1 on day
##1095.  Optimal control problem begins on *day* and calculates through 2 years post 
##introduction.  Used to examine the cost of delayed or preemptive vaccination programs.
##Results were not particularly interesting, but saving incase we want to revisit.
ebola.optim.timing<-function(params,M,day=0){
  
  # initial guesses - controls
  
  times<-seq(day,1825)
  inits=c(S1=params$N1,E1=0,I1=0,H1=0,D1=0,R1=0,
                 S2=params$N2,E2=0,I2=0,H2=0,D2=0,R2=0)
  
  v1 <- data.frame(times = times, v1 = rep(0,length(times)))
  v1_interp <- approxfun(v1, rule = 2)
  v1 = v1[,2]
  v2 = data.frame(times = times, v2 = rep(0,length(times)))
  v2_interp <- approxfun(v2, rule = 2)
  v2 = v2[,2]
  params<-c(params,v1_interp=v1_interp,v2_interp=v2_interp)
  # adjoints
  x = matrix(0, nrow = length(times), ncol = 13)
  lambda = matrix(0, nrow = length(times), ncol = 13)
  lambda_init = rep(0,12)
  names(lambda_init) = paste0("lambda",1:12)
  
  # define norm(X,1) command from matlab
  ##------------------------------------
  norm_oc <- function(x){colSums(abs(x))}
  
  delta<-.01
  counter <- 0
  test <- -1
  solx<-ode(y=inits,
            times=seq(1095,1825),
            func=ebola.ode,
            parms=params,
            method = "ode45")%>%
    data.frame()
  if(day<1095){
    solx<-bind_rows(data.frame(time=seq(day,1094),S1=params$N1,E1=0,I1=0,H1=0,D1=0,R1=0,
                               S2=params$N2,E2=0,I2=0,H2=0,D2=0,R2=0),
                    solx)
  }else if(day>1095){
    solx<-solx%>%filter(time>=day)
  }
  x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
  params<-c(params,x_interp=x_interp)
  while(test < 0 & counter < 200){
    counter <- counter + 1
    # set previous control, state, and adjoint 
    oldv1 <- v1
    oldv2 <- v2
    oldx <- solx
    oldlambda <- lambda
    
    # interpolate v
    params$v1_interp <- approxfun(times, v1, rule = 2)
    params$v2_interp <- approxfun(times, v2, rule = 2)
    
    
    # solve states
    if(day<1095){
      solx_1 <- ode(y = inits, times = seq(day,1095), func = ebola.ode, parms = params)
      solx_2 <- ode(y = solx_1[dim(solx_1)[1],-1]+c(-10,0,10,0,0,0,0,0,0,0,0,0), times = seq(1095,1825), func = ebola.ode, parms = params)
      solx<-bind_rows(data.frame(solx_1),data.frame(solx_2[-1,]))
      solx.forCosts<-solx
    }else if(day>=365){
      solx <- ode(y = inits+c(-10,0,10,0,0,0,0,0,0,0,0,0), times = seq(1095,1825), func = ebola.ode, parms = params)
      solx.forCosts<-solx
      solx<-solx%>%
        data.frame()%>%
        filter(time>=day)
      }
    
    
    
    # interpolate x (state)
    params$x_interp <- lapply(2:ncol(solx), function(x){approxfun(solx[,c(1,x)], rule = 2)})
    
    
    # solve adjoint
    lambda <- ode(y = lambda_init, times = rev(times), func = adjoints.ode, 
                  parms = params)
    
    lambda <- lambda[nrow(lambda):1,]
    
    # calculate v1* and v2*
    temp_v1 <- with(params,(-C1*(solx[,"S1"]+solx[,"I1"])+lambda[,"lambda1"]*solx[,"S1"]-lambda[,"lambda6"]*solx[,"S1"])/(2*epsilon1))
    temp_v2 <- with(params,(-C2*(solx[,"S2"]+solx[,"I2"])+lambda[,"lambda7"]*solx[,"S2"]-lambda[,"lambda12"]*solx[,"S2"])/(2*epsilon2))
    
    v1 <- pmin(M[1], pmax(0, temp_v1))
    v1 <- 0.5*(v1 + oldv1)
    v2 <- pmin(M[2], pmax(0, temp_v2))
    v2 <- 0.5*(v2 + oldv2)
    
    # recalculate test
    if(length(solx) != length(oldx)){browser()}
    
    test <- min(delta*norm_oc(data.frame(v1,v2))-norm_oc(data.frame(oldv1,oldv2)-data.frame(v1,v2)),
                delta*norm_oc(solx[,-1])-norm_oc(oldx[,-1]-solx[,-1]),
                delta*norm_oc(lambda[,-1])-norm_oc(oldlambda[,-1]-lambda[,-1]))
    print(counter)
    print(test)
    
  }
  dt<-.1
  J11=with(c(params,data.frame(solx.forCosts)),sum((b1*(betaI1*S1 + betaD1*S1*D1) )*dt))
  J12=with(c(params,data.frame(solx.forCosts)),sum((b2*(betaI2*S2*I2 + betaD2*S2*D2))*dt))
  J21=with(c(params,data.frame(solx.forCosts)),sum(C1*v1*(S1)+epsilon1*v1^2)*dt)
  J22=with(c(params,data.frame(solx.forCosts)),sum(C2*v2*S2+epsilon2*v2^2)*dt)
  solx<-data.frame(solx,v1,v2)
  returned<-list(solx,J11,J12,J21,J22)
  names(returned)<-c("SolX","J11","J12","J21","J22")
  return(returned)
}
