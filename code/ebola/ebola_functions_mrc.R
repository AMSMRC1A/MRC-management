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
  state<- sapply(1:length(x_interp),function(i){x_interp[[i]](t)})
  names(state) <- colnames(solx)[-1]
  v1<-v1_interp(t)
  v2<-v2_interp(t)
  
  with(as.list(c(Y,p,state)),{
    N1<-S1+E1+I1+H1+R1
    N2<-S2+E2+I2+H2+R2
    
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
