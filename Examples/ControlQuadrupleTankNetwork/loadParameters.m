function Cte = loadParameters()
%% Description
% This function gerenated a struct with constants of the model dynamics.
% Output:   -Cte: struct with the necessary constants

%% A FAZER
% POR ARGUMENO E SO FAZR LOAD DO QUE É IMPORTANTE PARA ACA COISA

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

%% Define algorithm parameters 
Cte.d = 10;%[1; 15; 1; 10];
Cte.T = 30;%[30; 30; 30; 30];
Cte.iLQRIt = 50;
Cte.iLQReps = 1e-4;
Cte.dTlin = 10;

%% LQR weigts
Cte.AntiWU = 10;

%% LQR weigts
Cte.R = eye(2);
Cte.h = eye(2,4);
Cte.H = zeros(4,6);
Cte.H(1:2,1:2) = eye(2);
Cte.H(3:end,5:end) = eye(2);
Cte.Qt = zeros(4);
Cte.Qi = 0.05*eye(2);
Cte.Qt(1:2,1:2) = 20*eye(2);
Cte.Qt(3:4,3:4) = 0.05*eye(2);
Cte.AntiWU = 10;
Cte.E = [1 0;
         0 1;
         0 0;
         0 0;
         1 0;
         0 1];
end