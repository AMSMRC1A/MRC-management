T=20;
R0=15;

test = -1;

delta = 0.001;
M = 1000;
tvec=linspace(0,T,M+1)';

x=zeros(M+1,4);
lambda=zeros(M+1,4);
u=zeros(M+1,1);
R=zeros(M+1,1);

while(test < 0)
    
    oldu = u;
    oldx = x;
    oldlambda = lambda;
    
    solx = ode45(@(t,x) statelab7(t,x,tvec,u),tvec,[1000 100 50 1165]);
    x = deval(solx,tvec)';

    sollamb = ode45(@(t,lambda) adjointlab7(t,lambda,tvec,x,u),[T 0],[0 0 0 0]);
    lambda = deval(sollamb,tvec)';
    
    S=x(:,1);
    lambda1=lambda(:,1);
    
    
    temp=(S.*lambda1)./2;
    u1 = min(0.9,max(0,temp));
    u = 0.5*(u1 + oldu);
    
    test=min([delta*norm(u,1)-norm(oldu-u,1) delta*norm(x,1)-norm(oldx-x,1) delta*norm(lambda,1)-norm(oldlambda-lambda,1)]);

end

R=x(:,4)-(x(:,2)+x(:,3)+x(:,1));

           subplot(3,2,1);plot(tvec,x(:,1))
           subplot(3,2,1);xlabel('Time')
           subplot(3,2,1);ylabel('S')
           subplot(3,2,2);plot(tvec,x(:,2))
           subplot(3,2,2);xlabel('Time')
           subplot(3,2,2);ylabel('E')
           subplot(3,2,3);plot(tvec,x(:,3))
           subplot(3,2,3);xlabel('Time')
           subplot(3,2,3);ylabel('I')
           subplot(3,2,4);plot(tvec,R)
           subplot(3,2,4);xlabel('Time')
           subplot(3,2,4);ylabel('R')
           subplot(3,2,5);plot(tvec,x(:,4))
           subplot(3,2,5);xlabel('Time')
           subplot(3,2,5);ylabel('N')
           subplot(3,2,6);plot(tvec,u)
           subplot(3,2,6);xlabel('Time')
           subplot(3,2,6);ylabel('u')
           subplot(3,2,6);axis([0 T -0.1 1])
           