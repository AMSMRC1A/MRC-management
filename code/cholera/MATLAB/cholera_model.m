%------------------------------------------------------------------
% cholera_model.m
%
% This script runs a model of cholera with optimal control for the
% 2021 MRC Management project.
%
% Instructions: See READ ME for steps to run the model.
% Required Files: cholera_states.m
%                 cholera_adjoints.m
%                 cholera_costs.m
%-------------------------------------------------------------------
clear all
close all

% Set the parameters for the simulation.
%mu1 = 0;        mu2 = 0;                    % Natural Birth/death rate
%beta_I1 = 2.14E-5;    beta_I2 = 2.14E-5;    % Direct Transmission
%beta_W1 = 1.01E-5;    beta_W2 = 1.01E-5;    % Water Transmission
%m1 = 0.05;         m2 = 0.05;               % Movement of healthy
%n1 = 0;         n2 = 0;                     % Movement of Ill
%gamma1 = 0.25;    gamma2 = 0.25;            % Recovery Rate
%delta1 = 5E-4;    delta2 = 5E-4;            % Death due to disease
%xi1 = 7.56E-3;       xi2 = 7.56E-3;         % Shedding into water
%nu1 = xi1;       nu2 = xi2;                 % Natural decay rate (equal to xi_i)
%rho1 = 0.05;      rho2 = 0.05;              % Flow of water downstream

% Set the parameters for the simulation.
mu1 = 0;            mu2 = 0;
%mu1 = 1E-4;        mu2 = 1E-4;                    % Natural Birth/death rate
beta_I1 = 0;       beta_I2 = 0;
%beta_I1 = 2.64E-5;    beta_I2 = 2.64E-5;    % Direct Transmission
beta_W1 = 1.21E-4;    beta_W2 = 1.21E-4;    % Water Transmission
m1 = 0.025;         m2 = 0.025;               % Movement of healthy 0.05
n1 = 0;         n2 = 0;                     % Movement of Ill 0.02
gamma1 = 0.25;    gamma2 = 0.25;            % Recovery Rate
delta1 = 5E-4;    delta2 = 5E-4;            % Death due to disease
xi1 = 7.56E-3;       xi2 = 7.56E-3;         % Shedding into water
nu1 = xi1;       nu2 = xi2;                 % Natural decay rate (equal to xi_i)
rho1 = 0.025;      rho2 = 0.025;              % Flow of water downstream 0.005 gave interesting

%Set up the parameters in the optimal control
%b1 = 1;     b2 = 1;
%C1 = .125;    C2 = .125;
%epsilon1 = 10;   epsilon2 = 10;
%M1 = 0.015;     M2 = 0.015;                % Maximum vaccination rate

%Set up the parameters in the optimal control
b1 = 5;     b2 = 1;
C1 = 0.625;    C2 = 0.125;
epsilon1 = 5E5;   epsilon2 = 1E5;
M1 = 0.015;     M2 = 0.015;                % Maximum vaccination rate

%Set up initial conditions for the model
%S1_init = 10000 - 100;      S2_init = 10000 - 10;
%I1_init = 100;              I2_init = 10;
%R1_init = 0;                R2_init = 0;
%W1_init = 100;              W2_init = 10;

S1_init = 9900;      S2_init = 9990;
I1_init = 100;              I2_init = 10;
R1_init = 0;                R2_init = 0;
W1_init = 100;              W2_init = 10;

%Set up time for the simulation
tfinal = 200;       tstep = 0.01; %Final time and time step for ODES

% Store the parameters for quick reference in member functions.
params = [mu1 mu2 beta_I1 beta_I2 beta_W1 beta_W2 m1 m2 n1 n2 gamma1 gamma2 delta1 delta2 xi1 xi2 nu1 nu2 rho1 rho2];
oc_params = [b1 b2 C1 C2 epsilon1 epsilon2];
ICs = [S1_init I1_init R1_init W1_init S2_init I2_init R2_init W2_init];
tvec = 0:tstep:tfinal;

% Set up the memory for states, adjoints, and controls
x = zeros(length(tvec), 8);
L = zeros(length(tvec), 8);
v1 = zeros(length(tvec), 1);
v2 = zeros(length(tvec), 1);

% Run the Optimal Control Algorithm
test = -1; 
tol = 0.001; %tolerance in updating test
count = 0;
fprintf('Count  | Test  |  |  Cost \n');
fprintf('------------------------- \n');
while (test < 1E-8)
    % Memory Rearrange
    oldv1 = v1;
    oldv2 = v2;
    oldv = [oldv1; oldv2];
    oldx = x;
    oldL = L;
    
    %Solve state ODE forward in time
    solx = ode45(@(t,x) cholera_states(t, x, tvec, v1, v2, params), tvec, ICs);
    x = deval(solx, tvec)';
    
    %Solve adjoing ODE backward in time
    solL = ode45(@(t, L) cholera_adjoints(t, L, x, tvec, v1, v2, params, oc_params), [tfinal, 0], zeros(8, 1));
    L = deval(solL, tvec)';
    
    %Calculate optimal v1 and v2
    temp_v1 = ((L(:, 1)-L(:,3)).*x(:, 1) - C1*x(:, 1))./(2*epsilon1);
    temp_v2 = ((L(:, 5) - L(:, 7)).*x(:, 5) - C2*x(:, 5))./(2*epsilon2);
    %Incorporate the bounds on controls
    v1 = min(M1, max(0, temp_v1));
    v2 = min(M2, max(0, temp_v2));
    
    %Update Controls
    v1 = 0.5*(v1+oldv1);
    v2 = 0.5*(v2+oldv2);
    
    v = [v1; v2];
    
    %Update test (stopping criterion)
    test = min([tol*norm(v, 1) - norm(oldv-v, 1) tol*norm(x, 1) - norm(oldx-x, 1) tol*norm(L, 1) - norm(oldL - L, 1)]);
    count = count+1;
    [J1b,J1c,J1e,J2b,J2c,J2e]=cholera_cost(oc_params, params, x, v1, v2, tfinal, tstep);
    fprintf(' %2d    %10.8f  %10.8f \n', count, test, J1b+J1c+J1e+J2b+J2c+J2e);
    
end
fprintf('Count  | Test  |  |  Cost \n');
fprintf('----------------------------------------------------------- \n');

%% No Control
%---------------------------------------------------------
% In this section, we solve the state equations to see the
% outbreak with no controls implemented.
%---------------------------------------------------------

%Set controls to 0
v1nC = zeros(length(tvec), 1);
v2nC = zeros(length(tvec), 1);

%Solve state equation
solxnC = ode45(@(t,x) cholera_states(t, x, tvec, v1nC, v2nC, params), tvec, ICs);
    xnC = deval(solxnC, tvec)';

%Find Costs for the system
[J1bnC,J1cnC,J1enC,J2bnC,J2cnC,J2enC] = cholera_cost(oc_params, params, xnC, v1nC, v2nC, tfinal, tstep);

%Display cost components
disp(strrep(['No Control Objective Functional for Patch 1 Infection = ' sprintf(' %d,', J1bnC) ''], ')', ')'))

disp(strrep(['No Control Objective Functional for Patch 1 Vaccination = ' sprintf(' %d,', J1cnC) ''], ')', ')'))

disp(strrep(['No Control Objective Functional for Patch 1 Epsilon = ' sprintf(' %d,', J1enC) ''], ')', ')'))
 
disp(strrep(['No Control Objective Functional for Patch 2 Infection = ' sprintf(' %d,', J2bnC) ''], ')', ')'))

disp(strrep(['No Control Objective Functional for Patch 2 Vaccination = ' sprintf(' %d,', J2cnC) ''], ')', ')'))

disp(strrep(['No Control Objective Functional for Patch 2 Epsilon = ' sprintf(' %d,', J2enC) ''], ')', ')'))

disp(strrep(['No Control Objective Functional for Both Patches = ' sprintf(' %d,', J1bnC+J1cnC+J1enC+J2bnC+J2cnC+J2enC) ''], ')', ')'))

fprintf('----------------------------------------------------------- \n');
%% Maximum control
%---------------------------------------------
% Now we will implement the maximum control
% in each patch and solve the state equations.
%---------------------------------------------

%Set Vaccinations to Maximum Rate
v1mC = M1*ones(length(tvec), 1);
v2mC = M2*ones(length(tvec), 1);

%Solve the state equations with maximum vaccination rate
solxmC = ode45(@(t,x) cholera_states(t, x, tvec, v1mC, v2mC, params), tvec, ICs);
    xmC = deval(solxmC, tvec)';

%Calculate costs of the system
[J1bmC,J1cmC,J1emC,J2bmC,J2cmC,J2emC] = cholera_cost(oc_params, params, xmC, v1mC, v2mC, tfinal, tstep);

%Display the cost components from the system
disp(strrep(['Maximum Control Objective Functional for Patch 1 Infection = ' sprintf(' %d,', J1bmC) ''], ')', ')'))

disp(strrep(['Maximum Control Objective Functional for Patch 1 Vaccination = ' sprintf(' %d,', J1cmC) ''], ')', ')'))

disp(strrep(['Maximum Control Objective Functional for Patch 1 Epsilon = ' sprintf(' %d,', J1emC) ''], ')', ')'))
 
disp(strrep(['Maximum Control Objective Functional for Patch 2 Infection = ' sprintf(' %d,', J2bmC) ''], ')', ')'))

disp(strrep(['Maximum Control Objective Functional for Patch 2 Vaccination = ' sprintf(' %d,', J2cmC) ''], ')', ')'))

disp(strrep(['Maximum Control Objective Functional for Patch 2 Epsilon = ' sprintf(' %d,', J2emC) ''], ')', ')'))

disp(strrep(['Maximum Control Objective Functional for Both Patches = ' sprintf(' %d,', J1bmC+J1cmC+J1emC+J2bmC+J2cmC+J2emC) ''], ')', ')'))

fprintf('----------------------------------------------------------- \n');
%% Optimal Control
%------------------------------------
% Display the components of the cost
% functional when the optimal control
% in each patch is used.
%------------------------------------
[J1b,J1c,J1e,J2b,J2c,J2e]=cholera_cost(oc_params, params, x, v1, v2, tfinal, tstep);

disp(strrep(['Objective Functional for Patch 1 Infection = ' sprintf(' %d,', J1b) ''], ')', ')'))

disp(strrep(['Objective Functional for Patch 1 Vaccination = ' sprintf(' %d,', J1c) ''], ')', ')'))

disp(strrep(['Objective Functional for Patch 1 Epsilon = ' sprintf(' %d,', J1e) ''], ')', ')'))
 
disp(strrep(['Objective Functional for Patch 2 Infection = ' sprintf(' %d,', J2b) ''], ')', ')'))

disp(strrep(['Objective Functional for Patch 2 Vaccination = ' sprintf(' %d,', J2c) ''], ')', ')'))

disp(strrep(['Objective Functional for Patch 2 Epsilon = ' sprintf(' %d,', J2e) ''], ')', ')'))

disp(strrep(['Objective Functional for Both Patches = ' sprintf(' %d,', J1b+J1c+J1e+J2b+J2c+J2e) ''], ')', ')'))

fprintf('----------------------------------------------------------- \n');

%% Uniform Approach
%-----------------------------------------------
% Now we implement a "one size fits all"
% approach in both patches. This is equivalent
% to implementing optimal control on the system
% where v1 = v2. In that case, we need to adjust
% the formulation of v
%-----------------------------------------------
% Set up the memory for states, adjoints, and controls
xu = zeros(length(tvec), 8);
Lu = zeros(length(tvec), 8);
vu = zeros(length(tvec), 1);
% Run the Optimal Control Algorithm
test = -1; 
tol = 0.001; %tolerance in updating test
count = 0;
fprintf('Count  | Test  |  |  Cost \n');
fprintf('------------------------- \n');
while (test < 1E-8)
    % Memory Rearrange
    oldvu = vu;
    oldxu = xu;
    oldLu = Lu;
    
    %Solve state ODE forward in time
    solxu = ode45(@(t,x) cholera_states(t, x, tvec, vu, vu, params), tvec, ICs);
    xu = deval(solxu, tvec)';
    
    %Solve adjoing ODE backward in time
    solLu = ode45(@(t, L) cholera_adjoints(t, L, x, tvec, vu, vu, params, oc_params), [tfinal, 0], zeros(8, 1));
    Lu = deval(solLu, tvec)';
    
    %Calculate optimal v1 and v2
    temp_vu = ((L(:, 1)-L(:,3)).*x(:, 1) - C1*x(:, 1)+(L(:, 5) - L(:, 7)).*x(:, 5) - C2*x(:, 5))./(2*(epsilon1+epsilon2));
    %Incorporate the bounds on controls
    vu = min(M1, max(0, temp_vu));
    
    %Update Controls
    vu = 0.5*(vu+oldvu);
    
    %Update test (stopping criterion)
    test = min([tol*norm(vu, 1) - norm(oldvu-vu, 1) tol*norm(xu, 1) - norm(oldxu-xu, 1) tol*norm(Lu, 1) - norm(oldLu - Lu, 1)]);
    count = count+1;
    [J1bu,J1cu,J1eu,J2bu,J2cu,J2eu]=cholera_cost(oc_params, params, x, vu, vu, tfinal, tstep);
    fprintf(' %2d    %10.8f  %10.8f \n', count, test, J1bu+J1cu+J1eu+J2bu+J2cu+J2eu);
    
end
fprintf('Count  | Test  |  |  Cost \n');
fprintf('----------------------------------------------------------- \n');

[J1bu,J1cu,J1eu,J2bu,J2cu,J2eu]=cholera_cost(oc_params, params, x, vu, vu, tfinal, tstep);

disp(strrep(['Uniform Functional for Patch 1 Infection = ' sprintf(' %d,', J1bu) ''], ')', ')'))

disp(strrep(['Uniform Functional for Patch 1 Vaccination = ' sprintf(' %d,', J1cu) ''], ')', ')'))

disp(strrep(['Uniform Functional for Patch 1 Epsilon = ' sprintf(' %d,', J1eu) ''], ')', ')'))
 
disp(strrep(['Uniform Functional for Patch 2 Infection = ' sprintf(' %d,', J2bu) ''], ')', ')'))

disp(strrep(['Uniform Functional for Patch 2 Vaccination = ' sprintf(' %d,', J2cu) ''], ')', ')'))

disp(strrep(['Uniform Functional for Patch 2 Epsilon = ' sprintf(' %d,', J2eu) ''], ')', ')'))

disp(strrep(['Uniform Functional for Both Patches = ' sprintf(' %d,', J1bu+J1cu+J1eu+J2bu+J2cu+J2eu) ''], ')', ')'))

fprintf('----------------------------------------------------------- \n');
%% %%%%%%%%%%%%%%%%%%%%    Plots    %%%%%%%%%%%%%%%%%
% This first set of plots displays all compartments for
% both patches and the case of no control, max control,
% or optimal control in each patch
%------------------------------------------------------

           subplot(2,4,1);
           hold on
           plot(tvec,x(:,1),'b-','LineWidth',2.5);
           plot(tvec,x(:,5),'r--','LineWidth',2.5)
           plot(tvec,xnC(:,1),'g-','LineWidth',1);
           plot(tvec,xnC(:,5),'g--','LineWidth',1)
           plot(tvec,xmC(:,1),'k-','LineWidth',1);
           plot(tvec,xmC(:,5),'k--','LineWidth',1)
           subplot(2,4,1);xlabel('Time')
           subplot(2,4,1);ylabel('S')
           
           subplot(2,4,2);
           hold on
           plot(tvec,x(:,2),'b-','LineWidth',2.5);
           plot(tvec,x(:,6),'r--','LineWidth',2.5)
           plot(tvec,xnC(:,2),'g-','LineWidth',1);
           plot(tvec,xnC(:,6),'g--','LineWidth',1)
           plot(tvec,xmC(:,2),'k-','LineWidth',1);
           plot(tvec,xmC(:,6),'k--','LineWidth',1)
           subplot(2,4,2);xlabel('Time')
           subplot(2,4,2);ylabel('I')
           
           subplot(2,4,3);
           hold on
           plot(tvec,x(:,3),'b-','LineWidth',2.5);
           plot(tvec,x(:,7),'r--','LineWidth',2.5)
           plot(tvec,xnC(:,3),'g-','LineWidth',1);
           plot(tvec,xnC(:,7),'g--','LineWidth',1)
           plot(tvec,xmC(:,3),'k-','LineWidth',1);
           plot(tvec,xmC(:,7),'k--','LineWidth',1)
           subplot(2,4,3);xlabel('Time')
           subplot(2,4,3);ylabel('R')
           
           subplot(2,4,4); 
           hold on
           plot(tvec,x(:,4),'b-','LineWidth',2.5);
           plot(tvec,x(:,8),'r--','LineWidth',2.5)
           plot(tvec,xnC(:,4),'g-','LineWidth',1);
           plot(tvec,xnC(:,8),'g--','LineWidth',1)
           plot(tvec,xmC(:,4),'k-','LineWidth',1);
           plot(tvec,xmC(:,8),'k--','LineWidth',1)
           legend('Patch 1 Optimal', 'Patch 2 Optimal', 'Patch 1 No Control','Patch 2 No Control','Patch 1 Max Control','Patch 2 Max Control')
           subplot(2,4,4);xlabel('Time')
           subplot(2,4,4);ylabel('W')
           
           subplot(2,4,5);
           hold on
           plot(tvec,v1,'LineWidth',2.5)
           plot(tvec,v1nC,'g-','LineWidth',1);
           plot(tvec,v1mC,'k-','LineWidth',1);
           legend('Patch 1 Optimal', 'Patch 1 No Control','Patch 1 Max Control')
           subplot(2,4,5);xlabel('Time')
           subplot(2,4,5);ylabel('v_1')
           subplot(2,4,5);axis([0 tfinal -0.0001 M1+.0001])  
              
           subplot(2,4,6);
           hold on
           plot(tvec,v2,'LineWidth',2.5)
           plot(tvec,v2nC,'g-','LineWidth',1);
           plot(tvec,v1mC,'k-','LineWidth',1);
           legend('Patch 2 Optimal', 'Patch 2 No Control','Patch 2 Max Control')
           subplot(2,4,6);xlabel('Time')
           subplot(2,4,6);ylabel('v_2')
           subplot(2,4,6);axis([0 tfinal -0.0001 M2+.0001])  

%% %%%%%%%%%%%%%%%%%%%%    Plots    %%%%%%%%%%%%%%%%%
% This second set of plots displays all compartments for
% both patches and the case of no control, max control,
% or optimal control in each patch
%------------------------------------------------------
figure;
           subplot(2,4,1);
           hold on
           plot(tvec,x(:,1),'r-','LineWidth',2.5);
           plot(tvec,x(:,5),'b-','LineWidth',2.5)
           plot(tvec,xu(:,1),'r--','LineWidth',1);
           plot(tvec,xu(:,5),'b--','LineWidth',1)
           subplot(2,4,1);xlabel('Time')
           subplot(2,4,1);ylabel('S')
           
           subplot(2,4,2);
           hold on
           plot(tvec,x(:,2),'r-','LineWidth',2.5);
           plot(tvec,x(:,6),'b-','LineWidth',2.5)
           plot(tvec,xu(:,2),'r--','LineWidth',1);
           plot(tvec,xu(:,6),'b--','LineWidth',1)
           subplot(2,4,2);xlabel('Time')
           subplot(2,4,2);ylabel('I')
           
           subplot(2,4,3);
           hold on
           plot(tvec,x(:,3),'r-','LineWidth',2.5);
           plot(tvec,x(:,7),'b-','LineWidth',2.5)
           plot(tvec,xu(:,3),'r--','LineWidth',1);
           plot(tvec,xu(:,7),'b--','LineWidth',1)
           subplot(2,4,3);xlabel('Time')
           subplot(2,4,3);ylabel('R')
           
           subplot(2,4,4); 
           hold on
           plot(tvec,x(:,4),'r-','LineWidth',2.5);
           plot(tvec,x(:,8),'b-','LineWidth',2.5)
           plot(tvec,xu(:,4),'r--','LineWidth',1);
           plot(tvec,xu(:,8),'b--','LineWidth',1)
           legend('Patch 1 Optimal', 'Patch 2 Optimal', 'Patch 1 Uniform','Patch 2 Uniform')
           subplot(2,4,4);xlabel('Time')
           subplot(2,4,4);ylabel('W')
           
           subplot(2,4,5);
           hold on
           plot(tvec,v1,'r-','LineWidth',2.5);
           plot(tvec,v2,'b-', 'LineWidth',2.5);
           plot(tvec,vu,'k--','LineWidth',1);
           legend('Patch 1 Optimal', 'Patch 2 Optimal','Uniform Control')
           subplot(2,4,5);xlabel('Time')
           subplot(2,4,5);ylabel('v')
           subplot(2,4,5);axis([0 tfinal -0.0001 M1+.0001])  
              