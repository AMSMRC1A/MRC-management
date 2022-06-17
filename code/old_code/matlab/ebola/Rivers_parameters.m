%%%%%%%%%%%%%%%%%%%%
%% Parameter file %%
%%%%%%%%%%%%%%%%%%%%


%% Initial conditions

% patch 1
N1 = 6e6;
S1 = N1 - 1;
E1 = 0;
I1 = 1;
H1 = 0;
D1 = 0;
R1 = 0;
ci = 0;

y0 = [S1, E1, I1, H1, D1, R1, ci];

%% Parameters
p.betaI = 0.16;     % transmission from I
p.betaH = 0.062;    % transmission from hospital
p.betaF = 0.489;    % transmission from funeral
p.alpha = 1/12;     % incubation period
p.gammaH = 1/3.24;  % time until hospitalization
p.gammaDH = 1/10;   % time from hospitalization to death
p.gammaF = 1/2;     % duration of traditional funeral
p.gammaI = 1/15;    % duration of infection
p.gammaD = 1/13;    % time from infection to death
p.gammaIH = 1/15;   % time from hospitalization to recovery
p.theta1 = 0.197;   % fraction of infected hospitalized
p.delta1 = 0.5;     % case fatality of unhospitalized
p.delta2 = 0.5;     % case fatality of hospitalized