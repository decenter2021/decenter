function Cte = getConstantsWL()
%% Description
% This function gerenated a struct with constants of the model dynamics.
% Output:   -Cte: struct with the necessary constants

%% Define algorithm parameters 
Cte.d = [1; 15; 1; 10];
Cte.T = [30; 30; 30; 30];
Cte.iLQRIt = 50;
Cte.iLQReps = 1e-4;
Cte.dTlin = 10;

%% Simulation parameters
Cte.dT = 1;
%% Define constants of the model

% Size of the system
Cte.n = 4;
Cte.m = 2;

% Sectional area of each tank
Cte.A1 = 28; %cm^2
Cte.A2 = 32; %cm^2
Cte.A3 = 28; %cm^2
Cte.A4 = 32; %cm^2

% Area of the hole in the bottom of each tank
Cte.a1 = 0.071; %cm^2
Cte.a2 = 0.057; %cm^2
Cte.a3 = 0.040; %cm^2
Cte.a4 = 0.040; %cm^2

Cte.A = [28;32;28;32];
Cte.a = [0.071; 0.057;0.04;0.04];
Cte.k = [3.33;3.33];
Cte.gamma = [0.7;0.6];

% % Sensivity of the water level sensors
% Cte.kc = 0.5; %V/cm

% Acceleration due to gravity
Cte.g = 981; %cm/s^2

% Pump flow sensivity
Cte.k1 = 3.33; %cm^3/Vs
Cte.k2 = 3.33; %cm^3/Vs
Cte.uMax = 12;

% Define the fraction of flow to the first and seconds tanks in the
% bifurcation valve
Cte.g1 = 0.7;
Cte.g2 = 0.6;

%% Define matrices 
Cte.alpha = [0.5041 1.7689 -1.8886; 1.13423 0.3249 -1.2141];
Cte.beta = [1.88883 -1.01092; -0.944417 1.76912];
  
%% LQR weigts

Cte.R = eye(Cte.m);

Cte.h = eye(Cte.n/2,Cte.n);

Cte.H = zeros(Cte.n,3*Cte.n/2);
Cte.H = zeros(Cte.n,3*Cte.n/2);
Cte.H(1:Cte.n/2,1:Cte.n/2) = eye(Cte.n/2);
Cte.H(Cte.n/2+1:end,Cte.n+1:end) = eye(Cte.m);

Cte.Qt = zeros(Cte.n);
Cte.Qi = 0.05*eye(Cte.n/2);
Cte.Qt(1:Cte.n/2,1:Cte.n/2) = 20*eye(Cte.n/2);
Cte.Qt(1+Cte.n/2:Cte.n,1+Cte.n/2:Cte.n) = Cte.Qi;

Cte.AntiWU = [10;10;10;10];

Cte.E = [1 0;
         0 1;
         0 0;
         0 0;
         1 0;
         0 1];

end