clear all
close all
clc

parameter_settings

T=100;

M1 = 0.05;
M2 = 0.05;

test = -1;

delta = 0.01;
M = 1000;
tvec=linspace(0,T,M+1)';

x=zeros(M+1,12);
lambda=zeros(M+1,12);
v1=zeros(M+1,1);
v2=zeros(M+1,1);


count = 1;

while(test < 0 && count < 50)
    
    oldv1 = v1;
    oldv2 = v2;
    oldx = x;
    oldlambda = lambda;
    
    solx = ode45(@(t,x) states(t,x,tvec,v1,v2,par),tvec,y0);
    x = deval(solx,tvec)';

    sollamb = ode45(@(t,lambda) adjoints(t,lambda,tvec,x,v1,v2,par),[T 0],zeros(1,12));
    lambda = deval(sollamb,tvec)';
    
    S1=x(:,1);
    E1=x(:,2);
    S2=x(:,7);
    E2=x(:,8);
    
    lambda1=lambda(:,1);
    lambda6=lambda(:,6);
    lambda7=lambda(:,7);
    lambda12=lambda(:,12);
    
    temp1= (-C1.*(S1 + E1) + lambda1.*S1 - lambda6.*S1)./(2*epsilon1);
    temp2= (-C2.*(S2 + E2) + lambda7.*S2 - lambda12.*S2)./(2*epsilon2);
    
    v11 = min(M1,max(0,temp1));
    v1 = 0.5*(v11 + oldv1);
    
    v21 = min(M2,max(0,temp2));
    v2 = 0.5*(v21 + oldv2);
    
    test=min([delta*norm(v1,1)-norm(oldv1-v1,1) delta*norm(v2,1)-norm(oldv2-v2,1) delta*norm(x,1)-norm(oldx-x,1) delta*norm(lambda,1)-norm(oldlambda-lambda,1)])

    count = count + 1
end

figure(1)
subplot(4,2,1);plot(tvec,x(:,1))
subplot(4,2,1);xlabel('Time')
subplot(4,2,1);ylabel('S1')
subplot(4,2,2);plot(tvec,x(:,2))
subplot(4,2,2);xlabel('Time')
subplot(4,2,2);ylabel('E1')
subplot(4,2,3);plot(tvec,x(:,3))
subplot(4,2,3);xlabel('Time')
subplot(4,2,3);ylabel('I1')
subplot(4,2,4);plot(tvec,x(:,4))
subplot(4,2,4);xlabel('Time')
subplot(4,2,4);ylabel('H1')
subplot(4,2,5);plot(tvec,x(:,5))
subplot(4,2,5);xlabel('Time')
subplot(4,2,5);ylabel('D1')
subplot(4,2,6);plot(tvec,x(:,6))
subplot(4,2,6);xlabel('Time')
subplot(4,2,6);ylabel('R1')
subplot(4,2,7);plot(tvec,v1)
subplot(4,2,7);xlabel('Time')
subplot(4,2,7);ylabel('v1')
subplot(4,2,7);axis([0 T -0.1 M1+0.05])

figure(3)
subplot(4,2,1);plot(tvec,lambda(:,1))
subplot(4,2,1);xlabel('Time')
subplot(4,2,1);ylabel('S1')
subplot(4,2,2);plot(tvec,lambda(:,2))
subplot(4,2,2);xlabel('Time')
subplot(4,2,2);ylabel('E1')
subplot(4,2,3);plot(tvec,lambda(:,3))
subplot(4,2,3);xlabel('Time')
subplot(4,2,3);ylabel('I1')
subplot(4,2,4);plot(tvec,lambda(:,4))
subplot(4,2,4);xlabel('Time')
subplot(4,2,4);ylabel('H1')
subplot(4,2,5);plot(tvec,lambda(:,5))
subplot(4,2,5);xlabel('Time')
subplot(4,2,5);ylabel('D1')
subplot(4,2,6);plot(tvec,lambda(:,6))
subplot(4,2,6);xlabel('Time')
subplot(4,2,6);ylabel('R1')

figure(2)
subplot(4,2,1);plot(tvec,x(:,7))
subplot(4,2,1);xlabel('Time')
subplot(4,2,1);ylabel('S2')
subplot(4,2,2);plot(tvec,x(:,8))
subplot(4,2,2);xlabel('Time')
subplot(4,2,2);ylabel('E2')
subplot(4,2,3);plot(tvec,x(:,9))
subplot(4,2,3);xlabel('Time')
subplot(4,2,3);ylabel('I2')
subplot(4,2,4);plot(tvec,x(:,10))
subplot(4,2,4);xlabel('Time')
subplot(4,2,4);ylabel('H2')
subplot(4,2,5);plot(tvec,x(:,11))
subplot(4,2,5);xlabel('Time')
subplot(4,2,5);ylabel('D2')
subplot(4,2,6);plot(tvec,x(:,12))
subplot(4,2,6);xlabel('Time')
subplot(4,2,6);ylabel('R2')
subplot(4,2,7);plot(tvec,v2)
subplot(4,2,7);xlabel('Time')
subplot(4,2,7);ylabel('v2')
subplot(4,2,7);axis([0 T -0.1 M2+0.05])

figure(4)
subplot(4,2,1);plot(tvec,lambda(:,7))
subplot(4,2,1);xlabel('Time')
subplot(4,2,1);ylabel('S2')
subplot(4,2,2);plot(tvec,lambda(:,8))
subplot(4,2,2);xlabel('Time')
subplot(4,2,2);ylabel('E2')
subplot(4,2,3);plot(tvec,lambda(:,9))
subplot(4,2,3);xlabel('Time')
subplot(4,2,3);ylabel('I2')
subplot(4,2,4);plot(tvec,lambda(:,10))
subplot(4,2,4);xlabel('Time')
subplot(4,2,4);ylabel('H2')
subplot(4,2,5);plot(tvec,lambda(:,11))
subplot(4,2,5);xlabel('Time')
subplot(4,2,5);ylabel('D2')
subplot(4,2,6);plot(tvec,lambda(:,12))
subplot(4,2,6);xlabel('Time')
subplot(4,2,6);ylabel('R2')
          