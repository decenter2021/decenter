%% Description
% This script simulates the extention of the original non-linear quadruple.
% The control vector is computed using predefined law (since we only want 
% to study the estimation performance). In this script observers of the LTV
% system are also implemented. It includes a centralized observer and two 
% distributed observers, whose gains are computed using the One Step and 
% Finite Horizon Algorithms proposed in [1]

%% Descentralization scheme
% The distributed scheme presented assumes the control vector is determined
% in a centralized manner and  is available to the local observers via 
% communication. It consists of 4 local observers in each tank,
% estimating their own position

%% Set constants and parameters that describe the network and simulation
clear;
% --- Network parameters ---
cte = getConstantsNTankNetwork(1);
% --- Simulation options ---
% Select initial state vector in cm
x0 = 20*ones(cte.n,1);
% Select the time window for the simulation in seconds
tspan = [0 800*cte.dT]; % (s)
% Compute number of discrete time-steps
simIt = tspan(2)/cte.dT+1;
LinIt = tspan(2)/cte.dTlin+1;
% Select default seed for random number generation, for debbuging purposes
rng(7); 

%% Initialize observe
% Define initial error covariance matrix
cte.P0 = 10*eye(cte.n);
xhat0 = x0+transpose(mvnrnd(zeros(cte.n,1),cte.P0));
% Initialize vector of observations
xhatC = zeros(cte.n,simIt-1);
xhatOS = zeros(cte.n,simIt-1);
% Initialize temporary cell to hold the covariance matrices
PC = cell(simIt,1);
POS = cell(simIt,1);
% Variable to hold the previous predicted estimation error covariance
PpredC = zeros(cte.n);
PpredOS = zeros(cte.n);

%% Discrete simulation for centralized and OS
% Discrete time vector
t_disc = 0:cte.dT:tspan(2);
% Initialize discrete control vector
uControlDiscrete = zeros(cte.m,simIt);
% Initialize measurement vector
y = cell(1,simIt-1);
x_disc = zeros(cte.n,simIt-1);
% Initialize the discrete time varying system that corresponds to the
% linearized non-linear model 
system = cell(simIt-1,7);
systemC = cell(LinIt-1,7);
systemOS = cell(LinIt-1,7);

for i = 1:simIt-1
    % Define indices for the linearized system model to take into account
    j_ = floor((i-1)/cte.dTlin)+1;   
    if i == 1     
         % Compute the linearized dynamics for the first instant 
         systemC(i,:) = getDiscreteDynamics(xhat0,cte);
         systemOS(i,:) = getDiscreteDynamics(xhat0,cte);
         % Get control vector with centralized LQR
         uControlDiscrete(:,floor(t_disc(i)/cte.dT)+1) = getControl((i-1)*cte.dT,cte.m);
         % Simulate nonlinear dynamics
         system(i,:) = getDiscreteDynamics(x0,cte);
         nonLinSol = ode45(@(t,x) xdotContinuous(x,uControlDiscrete(:,min(floor(t/cte.dT)+1,floor(t_disc(i+1)/cte.dT))),cte),[0 t_disc(i+1)],x0);
         x_disc(:,i) = deval(nonLinSol, t_disc(i+1)) + transpose(mvnrnd(zeros(cte.n,1),system{i,3}));
        
         % Observer dynamics
         % Measure output (systemC is used beacuse C and R are constant hence equal for systemC and systemOS )
         y{i,1} = systemC{j_,2}*x_disc(:,i) + transpose(mvnrnd(zeros(cte.o,1),systemC{j_,4}));
         
         % ----- Centralized -----
         % Compute observer gain
         PpredC = systemC{j_,1}*cte.P0*systemC{j_,1}'+systemC{j_,3};
         [KC,PpredC,PC{i,1}] = kalmanCentralizedLTV(systemC(j_,:),PpredC);
         % Predict
         xhatC_ = systemC{j_,6}+systemC{j_,1}*(xhat0-systemC{j_,6})+...
              systemC{j_,5}*(uControlDiscrete(:,floor(t_disc(i)/cte.dT)+1)-systemC{j_,7});
         % Update
         xhatC(:,i) = xhatC_ + KC*(y{i,1}-systemC{j_,2}*xhatC_);
         
         % ----- One Step -----
         % Compute one step observer gain
         PpredOS = systemOS{j_,1}*cte.P0*systemOS{j_,1}'+systemOS{j_,3};
         [KOS,PpredOS,POS{i,1}] = kalmanOneStepLTV(systemOS(j_,1:4),cte.E,PpredOS);
         % Predict
         xhatOS_ = systemOS{j_,6}+systemOS{j_,1}*(xhat0-systemOS{j_,6})+...
              systemOS{j_,5}*(uControlDiscrete(:,floor(t_disc(i)/cte.dT)+1)-systemOS{j_,7});
         % Update
         xhatOS(:,i) = xhatOS_ + KOS*(y{i,1}-systemOS{j_,2}*xhatOS_);
 
    else
         % Recompute the linearized dynamics of the system with a given 
         % periodicity
         if rem(i-1,cte.dTlin) == 0 
            systemC(j_,:) = getDiscreteDynamics(xhatC(:,i-1),cte);
            systemOS(j_,:) = getDiscreteDynamics(xhatOS(:,i-1),cte);
         end
         % Get control vector with predefined control law
         uControlDiscrete(:,floor(t_disc(i)/cte.dT)+1) = getControl((i-1)*cte.dT,cte.m);
         % Simulate nonlinear dynamics 
         [system(i,:)] = getDiscreteDynamics(x_disc(:,i-1),cte);
         nonLinSol = ode45(@(t,x) xdotContinuous(x,uControlDiscrete(:,min(floor(t/cte.dT)+1,floor(t_disc(i+1)/cte.dT))),cte),[t_disc(i) t_disc(i+1)],x_disc(:,i-1));
         x_disc(:,i) = deval(nonLinSol, t_disc(i+1)) + transpose(mvnrnd(zeros(cte.o,1),system{i,3}));
         % Observer dynamics
         % Measure output (systemC is used beacuse C and R are constant hence equal for systemC and systemOS )
         y{i,1} = systemC{j_,2}*x_disc(:,i) + transpose(mvnrnd(zeros(cte.n,1),systemC{j_,4}));
         % ----- Centralized -----
         % Compute observer gain
         [KC,PpredC,PC{i,1}] = kalmanCentralizedLTV(systemC(j_,1:4),PpredC);
         xhatC_ = systemC{j_,6}+systemC{j_,1}*(xhatC(:,i-1)-systemC{j_,6})+...
              systemC{j_,5}*(uControlDiscrete(:,floor(t_disc(i)/cte.dT)+1)-systemC{j_,7});
         xhatC(:,i) = xhatC_ + KC*(y{i,1}-systemC{j_,2}*xhatC_);
        
         % ----- One Step -----
         % Compute one step observer gain
         [KOS,PpredOS,POS{i,1}] = kalmanOneStepLTV(systemOS(j_,1:4),cte.E,PpredOS);
         % Predict
         xhatOS_ = systemOS{j_,6}+systemOS{j_,1}*(xhatOS(:,i-1)-systemOS{j_,6})+...
              systemOS{j_,5}*(uControlDiscrete(:,floor(t_disc(i)/cte.dT)+1)-systemOS{j_,7});
         % Update
         xhatOS(:,i) = xhatOS_ + KOS*(y{i,1}-systemOS{j_,2}*xhatOS_);
    end
end
  
%% Simulate FH propagated system
% Initialize the discrete time varying system that corresponds to the
% linearized non-linear model 
systemFH = cell(simIt-1,7);
% Initialize vector of predictions
x_predictFH = zeros(cte.n,simIt-1);
% Initialize vector of observations
xhatFH = zeros(cte.n,simIt-1);
% Define algorithm parameters 
FH_w = 40;
FH_d = 0;
opts.verbose = false;
opts.epsl = 1e-6;
opts.maxOLIt = 100;
% Compute number of FH windows
nW = floor((simIt-1)/(FH_w-FH_d));
% Number of remaining instants that do not fill a window
endW = rem(simIt-1,FH_w-FH_d); 
% Initialize temporary cell to hold the covariance matrices and gains
KFH = cell(simIt-1,1);
PFH = cell(simIt-1,1);

for i = 1:floor(nW)
    if i == 1
        %% First window
        % Compute system dynamics for the window
        systemFH(i,:) = getDiscreteDynamics(xhat0,cte);
        x_predictFH(:,1) = systemFH{i,1}*(xhat0-systemFH{i,6})+...
            systemFH{i,5}*(getControl((i-1)*cte.dT,cte.m)-systemFH{i,7})+systemFH{i,6};
        for j = 2:FH_w
            if rem(j-1,cte.dTlin) == 0 
                systemFH(j,:) = getDiscreteDynamics(x_predictFH(:,j-1),cte);
            else
                systemFH(j,:) = systemFH(j-1,:);
            end
            x_predictFH(:,j) = systemFH{j,1}*(x_predictFH(:,j-1)-systemFH{j,6})+...
                systemFH{j,5}*(getControl((j-1)*cte.dT,cte.m)-systemFH{j,7})+systemFH{j,6};
        end
        % Compute gains
        P0_pred = systemFH{i,1}*cte.P0*systemFH{i,1}'+systemFH{i,3};
        [auxK, auxP] = kalmanFiniteHorizonLTV(systemFH(1:end,1:4),cte.E,FH_w,P0_pred,opts);
        KFH((i-1)*(FH_w-FH_d)+1:i*(FH_w-FH_d),1) = auxK(1:(FH_w-FH_d));
        PFH((i-1)*(FH_w-FH_d)+1:i*(FH_w-FH_d),1) = auxP(1:(FH_w-FH_d));
        % Simulate dynamics
        for j = 1:FH_w-FH_d
             if j == 1
                 % Predict
                 xhatFH_ = systemFH{j,6}+systemFH{j,1}*(xhat0-systemFH{j,6})+...
                      systemFH{j,5}*(uControlDiscrete(:,floor(t_disc(j)/cte.dT)+1)-systemFH{j,7});
                 % Update
                 xhatFH(:,j) = xhatFH_ + KFH{1,1}*(y{j,1}-systemFH{j,2}*xhatFH_);
             else
                 % Predict
                 xhatFH_ = systemFH{j,6}+systemFH{j,1}*(xhatFH(:,j-1)-systemFH{j,6})+...
                      systemFH{j,5}*(uControlDiscrete(:,floor(t_disc(j)/cte.dT)+1)-systemFH{j,7});
                 % Update
                 xhatFH(:,j) = xhatFH_ + KFH{j,1}*(y{j,1}-systemFH{j,2}*xhatFH_);
             end 
        end
    else
        %% Remaining windows window
        % Compute system dynamics for the window
        j_ = (i-1)*(FH_w-FH_d)+1;
        systemFH(j_,:) = getDiscreteDynamics(xhatFH(:,j_-1),cte);
        x_predictFH(:,j_) = systemFH{j_,1}*(xhatFH(:,j_-1)-systemFH{j_,6})+...
            systemFH{j_,5}*(getControl((j_-1)*cte.dT,cte.m)-systemFH{j_,7})+systemFH{j_,6};
        for j = j_+1:j_+FH_w-1
            if rem(j-j_,cte.dTlin) == 0 
                systemFH(j,:) = getDiscreteDynamics(x_predictFH(:,j-1),cte);
            else
                systemFH(j,:) = systemFH(j-1,:);
            end
            x_predictFH(:,j) = systemFH{j,1}*(x_predictFH(:,j-1)-systemFH{j,6})+...
                systemFH{j,5}*(getControl((j-1)*cte.dT,cte.m)-systemFH{j,7})+systemFH{j,6};
        end
        % Compute gains
        P0_pred = systemFH{(i-1)*(FH_w-FH_d),1}*PFH{(i-1)*(FH_w-FH_d),1}*systemFH{(i-1)*(FH_w-FH_d),1}'+systemFH{(i-1)*(FH_w-FH_d),3};
        [auxK, auxP] = kalmanFiniteHorizonLTV(systemFH((i-1)*(FH_w-FH_d)+1:end,1:4),cte.E,min(FH_w,size(systemFH((i-1)*(FH_w-FH_d)+1:end,1),1)),P0_pred,opts);
        KFH((i-1)*(FH_w-FH_d)+1:i*(FH_w-FH_d),1) = auxK(1:(FH_w-FH_d));
        PFH((i-1)*(FH_w-FH_d)+1:i*(FH_w-FH_d),1) = auxP(1:(FH_w-FH_d));
        % Simulate dynamics
        for j = j_:j_+FH_w-FH_d-1
             % Predict
             xhatFH_ = systemFH{j,6}+systemFH{j,1}*(xhatFH(:,j-1)-systemFH{j,6})+...
                  systemFH{j,5}*(uControlDiscrete(:,floor(t_disc(j)/cte.dT)+1)-systemFH{j,7});
             % Update
             xhatFH(:,j) = xhatFH_ + KFH{j,1}*(y{j,1}-systemFH{j,2}*xhatFH_); 
        end

    end
end
%% Remaining partial window of the simulation 
if endW > 0
  % Compute system dynamics for the window
    j_ = i*(FH_w-FH_d)+1;
    systemFH(j_,:) = getDiscreteDynamics(xhatFH(:,j_-1),cte);
    x_predictFH(:,j_) = systemFH{j_,1}*(xhatFH(:,j_-1)-systemFH{j_,6})+...
            systemFH{j_,5}*(getControl((j_-1)*cte.dT,cte.m)-systemFH{j_,7})+systemFH{j_,6};
    for j = j_+1:j_-1+endW
        if rem(j-j_,cte.dTlin) == 0 
            systemFH(j,:) = getDiscreteDynamics(x_predictFH(:,j-1),cte);
        else
            systemFH(j,:) = systemFH(j-1,:);
        end
        x_predictFH(:,j) = systemFH{j,1}*(x_predictFH(:,j-1)-systemFH{j,6})+...
                systemFH{j,5}*(getControl((j-1)*cte.dT,cte.m)-systemFH{j,7})+systemFH{j,6};
    end
    % Compute gains
    P0_pred = systemFH{(floor(nW))*(FH_w-FH_d),1}*PFH{(floor(nW))*(FH_w-FH_d),1}*systemFH{(floor(nW))*(FH_w-FH_d),1}'+systemFH{(floor(nW))*(FH_w-FH_d),3};
    [auxK, auxP] = kalmanFiniteHorizonLTV(systemFH((i)*(FH_w-FH_d)+1:end,1:4),cte.E,endW,P0_pred,opts);  
    KFH(floor(nW)*(FH_w-FH_d)+1:end,1) = auxK(1:endW);
    PFH(floor(nW)*(FH_w-FH_d)+1:end,1) = auxP(1:endW);
    % Simulate dynamics
    for j = j_:j_-1+endW
         % Predict
         xhatFH_ = systemFH{j,6}+systemFH{j,1}*(xhatFH(:,j-1)-systemFH{j,6})+...
              systemFH{j,5}*(uControlDiscrete(:,floor(t_disc(j)/cte.dT)+1)-systemFH{j,7});
         % Update
         xhatFH(:,j) = xhatFH_ + KFH{j,1}*(y{j,1}-systemFH{j,2}*xhatFH_); 
    end
end
  

%% Plots
figure;
hold on;
grid on;
set(gca,'FontSize',30);
ylabel('$h_{13}$ (cm)','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
plot(t_disc,[x0(13) x_disc(13,:)], 'LineWidth',3);
plot(t_disc,[xhat0(13) xhatC(13,:)],'LineWidth',3);
plot(t_disc,[xhat0(13) xhatOS(13,:)],'LineWidth',3);
plot(t_disc,[xhat0(13) xhatFH(13,:)],'LineWidth',3);
legend('Nonlinear (true state)', 'Centralized', 'One-step', 'Finite-horizon (propagated)');
hold off;

figure;
hold on;
grid on;
set(gca,'FontSize',30);
ylabel('$h_{31}$ (cm)','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
plot(t_disc,[x0(31) x_disc(31,:)],'LineWidth',3);
plot(t_disc,[xhat0(31) xhatC(31,:)],'LineWidth',3);
plot(t_disc,[xhat0(31) xhatOS(31,:)],'LineWidth',3);
plot(t_disc,[xhat0(31) xhatFH(31,:)],'LineWidth',3);
legend('Nonlinear (true state)', 'Centralized', 'One-step', 'Finite-horizon (propagated)');
hold off;


%% getControl - Description
% This function outputs the known control law 
% Inputs:   - t: time instant
%           - m: number of pumps
% Outputs:  - du: control action 
function du = getControl(t,m)
    du = ones(m,1);
    if t < 400
        du = 6*du;
    else
        du = 4*du;
    end
end

%% getConstantsNTankNetwork - Description
% This function outputs a struct of constants of the model dynamics.
% Input:    -flagEqMatrices: only if true compute equilibrium matrices
% Output:   -Cte: struct with the necessary constants and parameters
function cte = getConstantsNTankNetwork(flagEqMatrices)
% --- Define dimentions of the model ---
cte.N = 40;      % number of tanks
% --- Size of the system ---
% number of states of the LTV system (not considering integral states)
cte.n = cte.N;  
cte.o = cte.N;    % number of inputs (pumps)
cte.m = cte.N/2;
% --- Define sampling intervals ---
cte.dT = 1;             % discrete sampling time (s)
cte.dTlin = 10;         % number of discrete sampling intervals MPC window
% --- Define parameters of the model ---
% Sectional area of each tank
cte.A = zeros(cte.N,1);
for i = 1:cte.N
    if rem(i,2) == 0
        cte.A(i) = 32; %cm^2
    else
        cte.A(i) = 28; %cm^2
    end
end
% Sectional area of the hole of each tank
cte.a = 0.040*ones(cte.N,1);
for i = 1:cte.N/2
    if rem(i,2) == 0
        cte.a(i) = 0.057; %cm^2
    else
        cte.a(i) = 0.071; %cm^2
    end
end
% Sensivity of the water level sensors
cte.kc = 0.5; %V/cm
% Acceleration due to gravity
cte.g = 981; %cm/s^2
% Pump flow sensivity
cte.k = 3.33*ones(cte.N/2,1);
% Fraction of flow to the bottom tanks in the bifurcation valve
cte.gamma = zeros(cte.N/2,1);
for i = 1:cte.N/2
    if rem(i,2) == 0
        cte.gamma(i) = 0.6; %cm^2
    else
        cte.gamma(i) = 0.7; %cm^2
    end
end
% --- Matrices for computing the alpha beta matrices ---
if flagEqMatrices % if true compute equilibrium matrices
    [alpha,beta] = getEquilibriumMatrices();
    cte.alpha = alpha;
    cte.beta = beta;
end
% --- Filter synthesis parameters ---
cte.R = 1*eye(cte.o);
% Choose Q
load('./data/NTanksNetworkEstimationOSFH_Q.mat','Q');
cte.Q = Q;
% Sparsity parttern
cte.E = eye(cte.n);

end

%% getEquilibriumMatrices - Description
% This function computes matrices alpha and beta according to [1] for
% equilibrium level computation
% Output:  -alpha, beta: matrices to compute water level 
% WARNING: Uses symbolic toolbox
function [alpha,beta] = getEquilibriumMatrices()
    cte = getConstantsNTankNetwork(0);  % get parameters 
    x = sym('x',[1 cte.n],'real');  % vector of positive water levels 
    assume(x,'positive');
    u = sym('u',[1,cte.m],'real');  % vector of positive pump actuations 
    assume(u,'positive');
    % vector of equations for null derivative of the water level
    xdot = sym(zeros(cte.n,1));    
    for i = 1:cte.n/2
       j = i+cte.n/2;
       xdot(i) = -(cte.a(i)/cte.A(i))*sqrt(2*cte.g)*x(i)+(cte.a(j)/cte.A(i))...
           *sqrt(2*cte.g*x(j))+cte.gamma(i)*cte.k(i)*u(i)/cte.A(i) == 0;
    end
    i = i+1;
    j = cte.n/2;
    xdot(cte.n/2+1) = -(cte.a(i)/cte.A(i))*sqrt(2*cte.g*x(i))+...
           (1-cte.gamma(j))*cte.k(j)*u(j)/cte.A(i) == 0;
    for i = cte.n/2+2:cte.n
       j = i-cte.n/2-1; 
       xdot(i) = -(cte.a(i)/cte.A(i))*sqrt(2*cte.g*x(i))+...
           (1-cte.gamma(j))*cte.k(j)*u(j)/cte.A(i) == 0;
    end
    % Solve xdot = 0 for equilibrium solution
    % Supress warning concerning vality of soltions according to
    % assunptions made in the variables
    warning('off')
    sol = solve(xdot,[x(cte.n/2+1:end) u(1:end)]);
    warning('on')
    sol = struct2cell(sol);
    % Plugin respective cooefficients in the entries of alpha and beta
    alpha = zeros(cte.n/2,cte.n/2+nchoosek(cte.n/2,2));
    for i = cte.n/2+1:cte.n
        [alpha(i-cte.n/2,:),~] = coeffs(sol{i-cte.n/2,1},x(1:cte.n/2));
    end
    beta = zeros(cte.m,cte.n/2);
    for i = cte.n+1:3*cte.n/2
        [beta(i-cte.n,:),~] = coeffs(sol{i-cte.n/2,1},x(1:cte.n/2));
    end
end

%% getDiscreteDynamics - Description
% This function computes the linearized discrete model for a given state.
% The equilibrium state is computed around the level in the lower tanks.
% Input:    - x: state vector
%           - cte: struct of constants of the model dynamics 
% Output:   - linDynamics: 1x7 cell with matrices A,C,Q,R,B, equilibrium
% state vector, and equilibrium control vector, uEq (in this order).
% Note: The state vector in this function includes integral states
function linDynamics = getDiscreteDynamics(x,cte)
% Initialize 1x7 cell to store the model matrices and equilibrium vectors
linDynamics = cell(1,7);
% Compute the equilibrium levels corresponding to the current water level 
% of the lower tanks, according to [1]
[xEq,uEq] = computeEquilibriumLevels(x(1:cte.n/2),cte);
% Compute linearized entries of matrices A and B for a continuous process
[A_cont,B_cont] = Dxdot(xEq,cte);
% Discretize dynamics with sampling interval cte.dT
G = expm([A_cont B_cont; zeros(cte.m,cte.n) zeros(cte.m,cte.m)]*cte.dT);
A = G(1:cte.n,1:cte.n);
B = G(1:cte.n,cte.n+1:end);
% Compute matrix A
linDynamics{1,1} = A;
% Matrix C
linDynamics{1,2} = cte.kc*eye(cte.n);
% Matrix Q
G = expm(cte.dT*[-A_cont cte.Q ; zeros(cte.n) A_cont']);
linDynamics{1,3} = transpose(G(cte.n+1:end,cte.n+1:end))*G(1:cte.n,cte.n+1:end);
% Matrix R
linDynamics{1,4} = cte.R;
% Compute matrix B for the system dynamics 
linDynamics{1,5} = B;
% Equilibrium states
linDynamics{1,6} = xEq;
linDynamics{1,7} = uEq;
end

%% computeEquilibriumLevels - Description
% This function computes the equilibrium levels corresponding to a given
% reference to the lower tanks, according to [1]
% Input :   - ref: reference vector to the lower tanks
%           - cte: parameters of the network 
% Output:   - xEq: vector of equilibrium water levels
%           - uEq: vector of equilibrium pump actuations
function [xEq,uEq] = computeEquilibriumLevels(ref,cte)
    % vector of reference levels and respective square root combinations
    xauxalpha = zeros(cte.n/2+nchoosek(cte.n/2,2),1);
    count = 1;
    for i = 1:cte.n/2
        for j = i:cte.n/2
            xauxalpha(count) = sqrt(ref(i))*sqrt(ref(j)); 
            count = count+1;
        end
    end
    % vector of the square root of the reference levels 
    xauxbeta = zeros(cte.n/2,1);
    for i = 1:cte.n/2
        xauxbeta(i) = sqrt(ref(i)); 
    end
    % Compute equilibrium levels and actuation
    xEq = zeros(cte.n,1);
    xEq(1:cte.n/2) = ref;
    xEq(cte.n/2+1:end) = cte.alpha*xauxalpha;
    uEq = cte.beta*xauxbeta;
end

%% xdotContinuous - Description
% This function computes the derivative of the state vector using the non 
% linear dynamics of the model, for a given state and control vector.
% Input:    - x: state vector 
%           - u: actuation vector
% Output:   - xdot: derivative of the state vector
function xdot = xdotContinuous(x,u,cte)
% Initialize xdot vector
xdot = zeros(cte.n,1);
% Compute xdot acconding to the non linear model
% For thr lower tanks
for i = 1:cte.n/2
   j = i+cte.n/2;
   xdot(i) = -(cte.a(i)/cte.A(i))*sqrt(2*cte.g*x(i))+(cte.a(j)/cte.A(i))*sqrt(2*cte.g*x(j))+...
       cte.gamma(i)*cte.k(i)*u(i,1)/cte.A(i);
end
% The first upper tank is connected to the same pump as the last lower tank
i = i+1;
j = cte.n/2;
xdot(cte.n/2+1) = -(cte.a(i)/cte.A(i))*sqrt(2*cte.g*x(i))+...
       (1-cte.gamma(j))*cte.k(j)*u(j)/cte.A(i);
% For the upper tanks except for the first
for i = cte.n/2+2:cte.n
   j = i-cte.n/2-1; 
   xdot(i) = -(cte.a(i)/cte.A(i))*sqrt(2*cte.g*x(i))+(1-cte.gamma(j))*cte.k(j)*u(j)/cte.A(i);
end      
end

%% Dxdot - Description 
% This function computes matrices A and B of the continuous linearized 
% model given an equilibrium state.
% Input:    - x_eq: equilibrium state vector
%           - cte: struct with constants of the model dynamics
function [A,B] = Dxdot(x_eq,cte)
% Initialize matrices
A = zeros(cte.n,cte.n);
B = zeros(cte.n,cte.m);
% Compute time constants of each tank
T = zeros(cte.n,1);
for i = 1:cte.n
   T(i) = (cte.A(i)/cte.a(i))*sqrt(2*x_eq(i)/cte.g);
end
% --- Compute A ---
% Lower tanks
for i = 1:cte.n
    A(i,i) = -1/T(i);
end
% Upper tanks
for i = 1:cte.n/2
    j = i+cte.n/2;
    A(i,j) = cte.A(j)/(cte.A(i)*T(j));
end
% --- Compute B ---
% Lower tanks
for i = 1:cte.n/2
    B(i,i) = cte.gamma(i)*cte.k(i)/cte.A(i);
end
% The pump of the last lower tanks feeds the first upper tank
B(cte.n/2+1,cte.n/2) = (1-cte.gamma(cte.n/2))*cte.k(cte.n/2)/cte.A(cte.n/2+1);
% Upper tanks except for the first
for i = cte.n/2+2:cte.n
    j = i-1-cte.n/2;
    B(i,j) =  (1-cte.gamma(j))*cte.k(j)/cte.A(i);
end
end

%% References
% [1] Pedroso L, Batista P, Oliveira P, Silvestre C. Discrete-time distributed
% Kalman filter design for networks of interconnected systems with linear 
% time-varying dynamics. International Journal of Systems Science. 2021; 
% https://doi.org/10.1080/00207721.2021.2002461

