%% Comparison of regulation performance between centralized, one-step, and
% finite-horizon methods - Monte Carlo simulations
% Methods developed in [1]
%% Generate random LTI system
n = 10;
m = 6;
rng(1); % Pseudo-random seed for consistency
% Alternatively comment out rng() to generate a random system
% Do not forget to readjust the synthesys parameters of the methods
A = rand(n,n)
B = rand(n,m)
Q = rand(n,n);
Q = Q*Q'
R = rand(m,m);
R = R*R'
E = round(rand(m,n))

%% Syntesize centralized, one-step, and finite-horizon controllers
[KC,PC] = LQRCentralizedLTI(A,B,Q,R);
[KOS,POS] = LQROneStepLTI(A,B,Q,R,E);
opts.W = 20;
opts.maxOLIt = 10;
[KFH,PFH] = LQRFiniteHorizonLTI(A,B,Q,R,E,opts);

%% Run Monte Carlo simulations
% Generate initial estimation error covariance accroding to the assumption
% P_0 = alpha*I
alpha = 1;
P0 = alpha*eye(n);
% Set number of Monte Carlo simulations
NMCSim = 1e4;
% Set simulation span of each Monte Carlo simulation
SimIt = 100;
% Initialise error cell, where the first dimension indicates the number of
% the Monte Carlo simulation; the second the instant of the simulation; and
% the third the method used
x = cell(NMCSim,SimIt,3);
u = cell(NMCSim,SimIt-1,3);
% This loops can, alternatively, be performed in parallel using parfor
for method = 1:3    % Iterate through the methods
    switch method   % Select gain
        case 1
            K = KC;
        case 2
            K = KOS;
        case 3 
            K = KFH;
    end
    for i = 1:NMCSim    % Iterate through the Monte Carlo simulations
        for k = 1:SimIt
            if k ~= 1
                x{i,k,method} = A*x{i,k-1,method}+B*u{i,k-1,method};
            else
                x{i,1,method} = transpose(mvnrnd(zeros(n,1),P0));
            end
            u{i,k,method} = -K*x{i,k,method};
        end    
    end
end

%% Finite window cost function computation
% Estimate expected values of J for each method
% Initialise Monte Carlo covariance cell, where the rows contain the data 
% of each method and the columns indicate the instant in the simulation
J = zeros(3,NMCSim);
% Iterate though the methods
for method = 1:3
    for i = 1:NMCSim
        J(method,i) = x{i,end,method}'*Q*x{i,end,method};
        for k = 1:SimIt-1
            J(method,i) = J(method,i) + ...
                x{i,k,method}'*Q*x{i,k,method}+u{i,k,method}'*R*u{i,k,method};            
        end 
    end
end
J = mean(J,2);

%% Results and comparison with projected performance
fprintf("Method\t\tCent.\tOS\tFH\n");
fprintf("Projected\t%.2f\t%.2f\t%.2f\n",...
    alpha*trace(PC),alpha*trace(POS),alpha*trace(PFH));
fprintf("Monte-Carlo\t%.2f\t%.2f\t%.2f\n",J(1),J(2),J(3));

%% References
% [1] Viegas D, Batista P, Oliveira P, Silvestre C. Distributed controller design 
% and performance optimization for discrete-time linear systems. Optim Control 
% Appl Meth. 2020;1-18. https://doi.org/10.1002/oca.2669
