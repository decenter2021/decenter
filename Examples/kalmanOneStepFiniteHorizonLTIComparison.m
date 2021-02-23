%% Comparison of estimation performance between centralized, one-step, and
% finite-horizon methods - Monte Carlo simulations
%% Generate random LTI system
n = 10;
o = 6;
rng(1); % Pseudo-random seed for consistency
% Alternatively comment out rng() to generate a random system
% Do not forget to readjust the synthesys parameters of the methods
A = rand(n,n)
C = rand(o,n)
Q = rand(n,n);
Q = Q*Q'
R = rand(o,o);
R = R*R'
E = round(rand(n,o))

%% Syntesize centralized, one-step, and finite-horizon filter gains
[KC,PC] = kalmanCentralizedLTI(A,C,Q,R);
[KOS,POS] = kalmanOneStepLTI(A,C,Q,R,E);
opts.W = 50;
opts.maxOLIt = 10;
[KFH,PFH] = kalmanFiniteHorizonLTI(A,C,Q,R,E,opts);

%% Run Monte Carlo simulations
% Generate initial estimation error covariance 
P0 = rand(n,n);
P0 = P0*P0';
% Set number of Monte Carlo simulations
NMCSim = 1e3;
% Set simulation span of each Monte Carlo simulation
SimIt = 100;
% Initialise error cell, where the first dimension indicates the number of
% the Monte Carlo simulation; the second the instant of the simulation; and
% the third the method used
error = cell(NMCSim,SimIt,3);
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
        error{i,1,method} = transpose(mvnrnd(zeros(n,1),P0));
        for k = 2:SimIt
            error{i,k,method} = (eye(n)-K*C)*(A*error{i,k-1,method}+...
                        transpose(mvnrnd(zeros(n,1),Q)))-...
                        K*transpose(mvnrnd(zeros(o,1),R));
        end
    end
end

%% Covariance matrix computation
% Estimate covariance matrix for each instant for each methods based on the
% error vectors of all the Monte Carlo simulations
% Initialise Monte Carlo covariance cell, where the rows contain the data 
% of each method and the columns indicate the instant in the simulation
P_MC = cell(3,SimIt);
% Iterate though the methods
for method = 1:3
    for k = 1:SimIt
        % Compute covariance matrix
        P_MC{method,k} = cov(cell2mat(error(:,k,method)')');
    end
end
% Compute the trace of each covariance matrix
traceEvolution = zeros(SimIt,3);
for method = 1:3
    for k = 1:SimIt
        traceEvolution(method,k) = trace(P_MC{method,k}); 
    end
end

%% Plot results 
% Plot the Trace vs instant of the simulation for the Monte Carlo
% simulation results
figure;
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
p = zeros(3,1);
for method = 1:3
    p(method) = plot(0:SimIt-1,traceEvolution(method,:),'LineWidth',3);
end
% Plot dashed projected performance 
plot(0:SimIt-1,trace(PC)*ones(1,SimIt),'--','LineWidth',2,'Color',get(p(1),'Color'));
plot(0:SimIt-1,trace(POS)*ones(1,SimIt),'--','LineWidth',2,'Color',get(p(2),'Color'));
plot(0:SimIt-1,trace(PFH)*ones(1,SimIt),'--','LineWidth',2,'Color',get(p(3),'Color'));
set(gcf, 'Position', [100 100 900 550]);
ylabel('$\mathrm{tr}(\mathbf{P}_{MC}(k))$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend('Cent.', 'OS', 'FH','Location','Best');
hold off;