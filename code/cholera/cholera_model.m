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

% Set the parameters for the simulation.
mu1 = 0;        mu2 = 0;                    % Natural Birth/death rate
beta_I1 = 2.14E-5;    beta_I2 = 2.14E-5;    % Direct Transmission
beta_W1 = 1.01E-5;    beta_W2 = 1.01E-5;    % Water Transmission
m1 = 0.05;         m2 = 0.05;               % Movement of healthy
n1 = 0;         n2 = 0;                     % Movement of Ill
gamma1 = 0.25;    gamma2 = 0.25;            % Recovery Rate
delta1 = 5E-4;    delta2 = 5E-4;            % Death due to disease
xi1 = 7.56E-3;       xi2 = 7.56E-3;         % Shedding into water
nu1 = xi1;       nu2 = xi2;                 % Natural decay rate (equal to xi_i)
rho1 = 0.05;      rho2 = 0.05;              % Flow of water downstream

%Set up the parameters in the optimal control
b1 = 1;     b2 = 1;
C1 = .125;    C2 = .125;
epsilon1 = 10;   epsilon2 = 10;
M1 = 0.015;     M2 = 0.015;                % Maximum vaccination rate

%Set up initial conditions for the model
S1_init = 10000 - 100;      S2_init = 10000 - 10;
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
while (test < 1E-16)
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
    fprintf(' %2d    %10.8f  %10.8f \n', count, test, cholera_cost(oc_params, params, x, v1, v2, tfinal, tstep));
    
end
