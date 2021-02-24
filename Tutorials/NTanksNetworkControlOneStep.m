%% Description
% This script simulates the extention of the original non-linear quadruple 
% tank network to a network of N tanks. Although the system is nonlinear,
% it can be linearized about equilibrium points. The known LTV
% decentralized control techniques can be applied to the linearized system.
% The control vector is computed making use of an MPC-like scheme, using 
% the one-step method and comparing its performance with the centralized
% solution. In addition, to take into consideration the nonlinearity of the
% dynamics of this network, an iLQR scheme is used.

%% Initilalize workspace
clear;
clear LQROneStepLTV; % clear permanent variables in LQROneStepLTV

%% Set constants and parameters that describe the network and simulation
% --- Network parameters ---
cte = getConstantsNTankNetwork(1);
% --- Simulation options ---
% Select initial state vector in cm
x0 = 20*ones(cte.n,1);
% Select the time window for the simulation in seconds
tspan = [0 600*cte.dT]; % (s)
% Compute number of discrete time-steps
simIt = tspan(2)/cte.dT+1;

%% Generate reference to lower tanks
% Initialize the reference vector for the lower tanks as a cell array,
% i.e., ref{1,k} : [h_ref1; ...; h_ref_N/2] for time instant k 
ref = cell(1,simIt);
% Initialize the equilibrum levels and inputs corresponding to the
% reference reference chosen for the lower tanks for each time instant
% i.e., refEq{1,k} : x_bar and refEq{1,k} : u_bar  for time instant k, 
% according to the control scheme proposed in [1].
% Note that intergral states are included in x_bar.
% Also note that a reference is needed for the instant after 
refEq = cell(2,simIt);
% Gererate random step and sine references with the defined periodicity,
% which repeat every four lower tanks
for k = 1:simIt
    % initalize vector of reference water levels for the lower tanks at
    % time instant k
    ref{1,k} = zeros(cte.n/2,1); 
    for j = 1:cte.n/2 % select one ammong 4 different reference waves
       switch rem(j,4)
           % sine wave with:
           %    - period: 50s
           %    - amplitude: 10cm
           %    - mean level: 30cm
           case 0   
               ref{1,k}(j) = 30+10*cos((k-1)/50);
           % square wave with:
           %    - period: 200s
           %    - amplitude: 5cm
           %    - mean level: 25cm
           case 1  
               if rem(floor((k-1)/100),int32(2)) == 0
                    ref{1,k}(j) = 30;
               else
                    ref{1,k}(j) = 20;
               end
           % sine wave with:
           %    - period: 35s
           %    - amplitude: 10cm
           %    - mean level: 30cm
           case 2
               ref{1,k}(j) = 30+10*cos((k-1)/35);
           % square wave with:
           %    - period: 400s
           %    - amplitude: 5cm
           %    - mean level: 25cm
           case 3
               if rem(floor((k-1)/200),int32(2)) == 0
                    ref{1,k}(j) = 30;
               else
                    ref{1,k}(j) = 20;
               end
       end
    end
    [xEq,uEq] = computeEquilibriumLevels(ref{1,k},cte);    
    refEq{1,k} = [xEq;zeros(cte.n/2,1)];
    refEq{2,k} = uEq;
end
% Do not allow square wave switching in the last time instant
ref{1,simIt}(1:2:cte.n/2,1) = ref{1,simIt-1}(1:2:cte.n/2,1);
[xEq,uEq] = computeEquilibriumLevels(ref{1,simIt},cte);    
refEq{1,simIt} = [xEq;zeros(cte.n/2,1)];
refEq{2,simIt} = uEq;

%% Discrete simulation and control
% Discrete time vector
t_disc = 0:cte.dT:tspan(2);
% Initialize temporary cell to hold the gain matrices for each of the four
% methods used
K = cell(4,simIt-1);
% Initialize discrete control vector
uControlDisc = cell(4,simIt-1);
% Initialize cell array to hold the component of the command action u_a
u_a = cell(4,simIt-1);
% Initialize measurement and state vectors
xNonLin = cell(4,simIt);
x = cell(4,simIt);
% Initialize nonlinear solver for each control action
nonLinSol = cell(4,1);
% --- Brief details of the control laws simulated (see [1] for more) --
% An MPC sheme is used where, for each time-instant, a window of gains into
% the future is computed. Then, only a fraction of those gains are actually
% used in the pumps. The tracker design is put forward in [1]. As the
% conditions of [Remark 5, 1] are satisfied, that particular case of
% [Theorem 2, 1].
% Four distinct control laws are compared:
% 1. Centralized LQR gains (d=1 gains are used out of each MPC window)
% 2. Centralized LQR gains (d>1 gains are used out of each MPC window)
% d is defined in private function getConstantsNTankNetwork 
% 3. One-step LQR gains (d=1 gains are used out of each MPC window)
% 4. One-step LQR gains (d>1 gains are used out of each MPC window)
% d is defined in private function getConstantsNTankNetwork 
% ---
% Iterate through every discrete time-instant
for k = 1:simIt
    if k == 1   % For the first instant initilaization is necessary
        for m = 1:4 % Iterate through every method
            % --- Measure output ---
            % water level obtained in each water level sensor (wo/ noise)
            x{m,k}(1:cte.n,1) = x0;
            % Initialize intergrator state
            x{m,k}(cte.n+1:3*cte.n/2) = x0(1:cte.n/2,1)-ref{1,k}(1:cte.n/2,1);
            % Initialize nonliner solver 
            xNonLin{m,k} = x0;
            % Compute the first MPC window
            [K(m,k:min(k+cte.d(m)-1,simIt-1)),u_a(m,k:min(k+cte.d(m)-1,simIt-1))]...
                = iLQR(x{m,k},refEq(:,k:end),cte,min(cte.d(m),simIt-k),m);
            % Control law [1]: u(k) = -K(k)(x(k)-x_bar(k))+u_bar(k)+u_a(k)
            uControlDisc{m,k} = -K{m,k}*(x{m,k}-refEq{1,k})+u_a{m,k}+refEq{2,k};
            % Saturate commands to the pumps
            uControlDisc{m,k}(uControlDisc{m,k}< cte.uMin) = cte.uMin;
            uControlDisc{m,k}(uControlDisc{m,k}>cte.uMax) = cte.uMax;
            % Simulate nonlinear dynamics
            nonLinSol{m,1} = ode45(@(t,x) xdotContinuous(x,uControlDisc{m,min(floor(t/cte.dT)+1,...
                round(t_disc(k+1)/cte.dT))},cte),[t_disc(k) t_disc(k+1)],x0);
        end
    elseif k ~= simIt
         for m = 1:4 % Iterate through every method
             % --- Measure output ---
             x{m,k}(1:cte.n,1) = deval(nonLinSol{m,1}, t_disc(k));
             % --- Integrate tracking error ---
             x{m,k}(cte.n+1:3*cte.n/2) = x{m,k-1}(cte.n+1:3*cte.n/2)+x{m,k}(1:cte.n/2,1)...
                 -ref{1,k}(1:cte.n/2,1);
             % Anti windup for integral action according to [1]
             for j = 1:cte.n/2
                if abs(x{m,k}(cte.n+j)) > cte.AntiWU(m)
                    x{m,k}(cte.n+j) = cte.AntiWU(m)*abs(x{m,k}(cte.n+j))/x{m,k}(cte.n+j); 
                end
             end
             xNonLin{m,k} = deval(nonLinSol{m,1}, t_disc(k));
             % If a new MPC window needs to be computed this time-instant
             if rem(k-1,cte.d(m)) == 0
                [K(m,k:min(k+cte.d(m)-1,simIt-1)),u_a(m,k:min(k+cte.d(m)-1,simIt-1))]...
                    = iLQR(x{m,k},refEq(:,k:end),cte,min(cte.d(m),simIt-k),m);
             end

             uControlDisc{m,k} = -K{m,k}*(x{m,k}-refEq{1,k})+u_a{m,k}+refEq{2,k};
             uControlDisc{m,k}(uControlDisc{m,k}<0) = 0;
             uControlDisc{m,k}(uControlDisc{m,k}>cte.uMax) = cte.uMax;
             % Simulate nonlinear dynamics
             nonLinSol{m,1} = ode45(@(t,x) xdotContinuous(x,uControlDisc{m,min(floor(t/cte.dT)+1,round(t_disc(k+1)/cte.dT))},cte),[t_disc(k) t_disc(k+1)],deval(nonLinSol{m,1}, t_disc(k)));
         end
    else % for the last time-instant
        % compute only the state at the last time-instant
        for m = 1:4
            xNonLin{m,k} = deval(nonLinSol{m,1}, t_disc(k));
        end
    end
    % Show status
    if k == 1 
        fprintf('Running simulation of N tank network: '); 
        delstatus = ''; 
    end
    status = strcat(sprintf('%3.1f', 100*k/simIt),'%%');
    fprintf([delstatus, status]);
    delstatus= repmat(sprintf('\b'),1,length(status)-1);
    if k == simIt fprintf('\n'); end
end
%% Plot results
% Plot evolution of water levels
for tank = 1:cte.n
    figure;
    xAux = zeros(cte.n,simIt);
    refAux = zeros(1,simIt);
    for k = 1:simIt
        for m = 1:4
            xAux(m,k) = xNonLin{m,k}(tank,1);
            if tank <= cte.n/2
                refAux(1,k) = ref{1,k}(tank,1);
            end
        end
    end
    hold on;
    set(gca,'FontSize',35);
    ylabel(sprintf("$h_{%d}$ (cm)",tank),'Interpreter','latex');
    xlabel('$t$ (s)','Interpreter','latex');
    for m = 1:4
        p = plot(t_disc,xAux(m,:));
        p.LineWidth = 3;
    end
    if tank <= cte.n/2
         p = plot(t_disc,refAux(1,:));
         p.LineWidth = 3;
        legend(sprintf('Centralized (d = %d)',cte.d(1)),sprintf('Centralized (d = %d)',cte.d(2)),...
            sprintf('One-step (d = %d)',cte.d(3)),sprintf('One-step (d = %d)',cte.d(4)),'Reference');
    else
        legend(sprintf('Centralized (d = %d)',cte.d(1)),sprintf('Centralized (d = %d)',cte.d(2)),...
            sprintf('One-step (d = %d)',cte.d(3)),sprintf('One-step (d = %d)',cte.d(4)));
    end
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    hold off;   
end
% Plot evolution of command action
for pump = 1:cte.n/2
    figure;
    uAux = zeros(4,simIt-1);
    for k = 1:simIt-1
        for m = 1:4
            uAux(m,k) = uControlDisc{m,k}(pump,1);
        end
    end
    hold on;
    set(gca,'FontSize',35);
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ylabel(sprintf("$u_{%d}$ (V)",pump),'Interpreter','latex');
    xlabel('$t$ (s)','Interpreter','latex');
    for m = 1:4
        p = plot(t_disc(1:end-1),uAux(m,:));
        p.LineWidth = 3;
    end
    legend(sprintf('Centralized (d = %d)',cte.d(1)),sprintf('Centralized (d = %d)',cte.d(2)),...
            sprintf('One-step (d = %d)',cte.d(3)),sprintf('One-step (d = %d)',cte.d(4)));
end

%% iLQR - Description
% This function computes the iLQR (iterative LQR gains). It is necessary
% because of the nonlinearity of the N tank network. The computation of the
% gains for a finite window using either a centralized or decentralized
% method require that the future dynamics of the system are known. However,
% for this system, the system dynamics in the future depend on the future
% state. iLQR iterative procedure where the future dynamics are being
% updated every time a new MPC window is computed, until convergence. For
% more details see [1].
% Input:    - x0: state at the beginning of the new window
%           - ref: cell array containing x_bar and u_bar
%           - d: number of gains to output out of those computed for the
%           whole window
%           - cte: struct with the necessary constants and parameters
%           - m: number of the method used for the computation of LQR gains
% Output:   - K: struct of LQR gains
%           - u_a: struct of addictional command actions
function [K,u_a] = iLQR(x0,ref,cte,d,m)
    % Initialize cell array for the future dynamics of the system 
    system = cell(cte.T(m)+1,7);
    % --- Forward pass variables ---
    % sequence of states throughout the window
    x = cell(1,cte.T(m)+1);
    % sequence of additional command action throughout the window
    u_a = cell(1,cte.T(m)+1); 
    % sequence of command action throughout the window
    uControlDisc = cell(1,cte.T(m));
    % sequence of command action throughout the window (previous iteration)
    % to check when convergence is reached
    prevuControlDisc = cell(1,cte.T(m));
    % Time of discrete-time intants
    t_disc = 0:cte.dT:cte.T(m);
    % Perform the iLQR iterations up to a maximum of cte.iLQRIt iterations
    for k = 1:cte.iLQRIt
       % Forward pass
       if k == 1
           % Initial propagated system assumes level is mantained constant
           system(1,:) = getDiscreteDynamics(x0,cte);
           for i = 2:cte.T(m)+1
               system(i,:) = system(1,:);
           end
       else
           for i = 1:cte.T(m) % Simulate forward pass
                if i == 1
                    % Compute the linearized dynamics for the first instant 
                    system(i,:) = getDiscreteDynamics(x0,cte);
                    % (assumes level is mantained constant)
                    for l = i+1:i+cte.dTlin-1
                        system(l,:) = system(i,:);
                    end
                    % Measure simulated output
                    x{1,i} = x0;
                    % Compute aditional command action
                    u_a{1,i} = (cte.h*system{i,2}(1:cte.n,1:cte.m))\cte.h*...
                       (ref{1,min(i+1,size(ref,2))}(1:cte.n,1)-ref{1,min(i,size(ref,2))}(1:cte.n,1));
                    % Control law [1]: u(k) = -K(k)(x(k)-x_bar(k))+u_bar(k)+u_a(k)
                    uControlDisc{1,i} = -K{i,1}*(x{1,i}-ref{1,min(i,size(ref,2))})...
                        +u_a{1,i}+ref{2,min(i,size(ref,2))};
                    % Saturate commands to the pumps
                    uControlDisc{1,i}(uControlDisc{1,i}<cte.uMin) = cte.uMin;
                    uControlDisc{1,i}(uControlDisc{1,i}>cte.uMax) = cte.uMax;
                    % Simulate nonlinear dynamics
                    nonLinSol = ode45(@(t,x) xdotContinuous(x,uControlDisc{1,min(floor(t/cte.dT)+1,round(t_disc(i+1)/cte.dT))},cte),[t_disc(i) t_disc(i+1)],x0(1:cte.n,1));
                    x{1,i+1}(1:cte.n,1) = deval(nonLinSol, t_disc(i+1));
                    % Integrate tracking error
                    x{1,i+1}(cte.n+1:3*cte.n/2) = x{1,i}(cte.n+1:3*cte.n/2)+x{1,i+1}(1:cte.n/2,1)-ref{1,i+1}(1:cte.n/2,1);
                else
                    % if a new linearization is necessary
                    if rem(i-1,cte.dTlin) == 0 
                       system(i,:) = getDiscreteDynamics(x{1,i},cte);
                       for l = i+1:cte.T(m)+1
                           system(l,:) = system(i,:);
                       end
                    end
                    % Compute aditional command action
                    u_a{1,i} = (cte.h*system{i,2}(1:cte.n,1:cte.m))\cte.h*(ref{1,min(i+1,size(ref,2))}(1:cte.n,1)-ref{1,min(i,size(ref,2))}(1:cte.n,1));
                    % Control law [1]: u(k) = -K(k)(x(k)-x_bar(k))+u_bar(k)+u_a(k)
                    uControlDisc{1,i} = -K{i,1}*(x{1,i}-ref{1,min(i,size(ref,2))})+u_a{1,i}++ref{2,min(i,size(ref,2))};
                    uControlDisc{1,i}(uControlDisc{1,i}<cte.uMin) = cte.uMin;
                    uControlDisc{1,i}(uControlDisc{1,i}>cte.uMax) = cte.uMax;
                    % Simulate nonlinear dynamics
                    nonLinSol = ode45(@(t,x) xdotContinuous(x,uControlDisc{1,min(floor(t/cte.dT)+1,round(t_disc(i+1)/cte.dT))},cte),[t_disc(i) t_disc(i+1)],x{1,i}(1:cte.n,1));
                    x{1,i+1}(1:cte.n,1) = deval(nonLinSol, t_disc(i+1));
                    % Integrate tracking error
                    x{1,i+1}(cte.n+1:3*cte.n/2) = x{1,i}(cte.n+1:3*cte.n/2)+x{1,i+1}(1:cte.n/2,1)-ref{1,min(i+1,size(ref,2))}(1:cte.n/2,1);
                    % Anti windup for integral action according to [1]
                    for j = 1:cte.n/2
                       if abs(x{1,i}(cte.n+j)) > cte.AntiWU(m)
                           x{1,i}(cte.n+j) = cte.AntiWU(m)*abs(x{1,i}(cte.n+j))/x{1,i}(cte.n+j); 
                       end
                    end
                end
           end
       end
       % --- stopping criterion --- 
       % stop the iteartions if the maximum difference
       % in relation to the actuation computed in the previous iteration
       % falls under cte.iLQReps
       if k > 2
           dif = zeros(1,cte.T(m));
           for i = 1:cte.T(m)
                dif(1,i) = norm(prevuControlDisc{1,i}-uControlDisc{1,i})/...
                    norm(uControlDisc{1,i});
           end
           if max(dif) < cte.iLQReps 
               break; 
           end
       end
       prevuControlDisc = uControlDisc;
       % --- Compute LQR gains ---
       if m<= 2 % Centralized
            [K,~] = LQRCentralizedLTV(system(:,1:4),cte.T(m));
       else % One-step
            [K,~] = LQROneStepLTV(system(:,1:4),cte.E,cte.T(m));
       end
       % check if maximum number of iterations was reached and issue
       % warning
       if k == cte.iLQRIt
           fprintf(sprintf('The maximum number of iLQR iterations was reached before convergence for method number %d.\n',m));
       end
    end
    % Output only gains and additional command action that are used
    K = transpose(K(1:d,1));
    u_a = u_a(1,1:d);
end

%% getConstantsNTankNetwork - Description
% This function outputs a struct of constants of the model dynamics.
% Input:    -flagEqMatrices: only if true compute equilibrium matrices
% Output:   -Cte: struct with the necessary constants and parameters
function cte = getConstantsNTankNetwork(flagEqMatrices)
% --- Define dimentions of the model ---
cte.N = 4;      % number of tanks
% --- Size of the system ---
% number of states of the LTV system (not considering integral states)
cte.n = cte.N;  
cte.m = cte.N/2;    % number of inputs (pumps)
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
cte.uMax = 12;
cte.uMin = 0;
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
% --- LQR synthesis parameters ---
% Command action weighting matrix
cte.R = eye(cte.m);
% Output tracking matrix
cte.h = eye(cte.n/2,cte.n); % output tracking matrix wo/ integral action
cte.H = zeros(cte.n,3*cte.n/2);
cte.H = zeros(cte.n,3*cte.n/2);
cte.H(1:cte.n/2,1:cte.n/2) = eye(cte.n/2);
cte.H(cte.n/2+1:end,cte.n+1:end) = eye(cte.m);
% Tracking performance and intergral action weighting matrices
cte.Qt = zeros(cte.n);
cte.Qi = 0.05*eye(cte.n/2);
cte.Qt(1:cte.n/2,1:cte.n/2) = 20*eye(cte.n/2);
cte.Qt(1+cte.n/2:cte.n,1+cte.n/2:cte.n) = cte.Qi;
% Anti wind-up for integrator
cte.AntiWU = [10;10;10;10];
% Sparsity parttern
cte.E = zeros(cte.m,3*cte.n/2);
cte.E(1:cte.m,1:cte.n/2) = eye(cte.n/2);
cte.E(1:cte.m,cte.n+1:end) = eye(cte.n/2);

% --- Define algorithm parameters ---
% (vector of parameters for each method), i.e., cte.d(1) = d_method1,...
cte.d = [1; 15; 1; 10];     % number of finite-window gains used
cte.T = [30; 30; 30; 30];   % finite-window length
cte.iLQRIt = 100;           % maximum number of iterations of the iLQR alg
cte.iLQReps = 1e-4;         % parameter for the stopping criterion of iLQR
end

%% getConstantsNTankNetwork - Description
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
[xEq,uEq] = computeEquilibriumLevels(x(1:cte.n/2,1),cte);
% Compute linearized entries of matrices A and B for a continuous process
[A_cont,B_cont] = Dxdot(xEq,cte);
% Discretize dynamics with sampling interval cte.dT
G = expm([A_cont B_cont; zeros(cte.m,cte.n) zeros(cte.m,cte.m)]*cte.dT);
A = G(1:cte.n,1:cte.n);
B = G(1:cte.n,cte.n+1:end);
% Compute matrix A for the system dynamics including integral action, 
% according to [1]
linDynamics{1,1}(1:cte.n,1:cte.n) = A;
linDynamics{1,1}(cte.n+1:3*cte.n/2,1:cte.n) = cte.h*A;
linDynamics{1,1}(1:cte.n,cte.n+1:3*cte.n/2) = zeros(cte.n,cte.n/2);
linDynamics{1,1}(cte.n+1:3*cte.n/2,cte.n+1:3*cte.n/2) = eye(cte.n/2);
% Compute matrix B for the system dynamics including integral action, 
% according to [1]
linDynamics{1,2}(1:cte.n,1:cte.m) = B;
linDynamics{1,2}(cte.n+1:3*cte.n/2,1:cte.m) = cte.h*B;
% State weigting matrix (icncluding integral action)
linDynamics{1,3} = cte.H'*cte.Qt*cte.H;
% Command action weighting matrix
linDynamics{1,4} = cte.R;
% Equilibrium states
linDynamics{1,5} = [xEq;zeros(cte.n/2,1)];
linDynamics{1,6} = uEq;
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
% [1] L. Pedroso, and P. Batista (xxx), Discrete-time decentralized linear 
% quadratic control for linear time-varying systems, Int J Robust Nonlinear
% Control, xxx;xx:x?x. [Submitted to journal]
