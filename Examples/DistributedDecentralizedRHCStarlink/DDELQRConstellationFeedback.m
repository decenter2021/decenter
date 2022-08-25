%% Simulation of constellation of LEO satellites
% Implementation of distributed decentralized LQR

%% Notes and issues:
% - 'Define constellation' section should be uploaded from .h file when it
% is implemented with tudat
% - allow a starting time different from 0 / and an associated anchor to
% allow for the continuation of a simulation

%% Init 
% Variable 'matlab_check' may be received as an input to run a small
% simulation propagated in matlab to make sure there are no errors
clearvars -except matlab_check; % Clear workspace variables
% Add source of orbital dynamics utilities
addpath('./src');
addpath('../src-osculating2mean');
% Add source of tudat feedback matlab server class
addpath('../src-tudat'); 
tic; % Log the execution time

%% Define execution options 
% Turn off tudat simulation for debug
tudatSimulation = true;
% Turn off parameter upload from C header
uploadHeaderParameters = true;
headerParametersFilepath = 'constellationParameters.h';
% Turn off decentralized computation for debug purposes
decentralized = true;
% If the goal is to perform a matlab check propagate dynamics in matlab and
% upload paramenters for C header file
if exist('matlab_check','var') && matlab_check
    tudatSimulation = false;
    uploadHeaderParameters = true;
end

%% Initialize server
if tudatSimulation
    if ~uploadHeaderParameters
        error('@MATLAB server: Parameters must be uploaded from C header file for tudat simulation.');
    end
fprintf('@MATLAB server: Setting up feedback server.\n');
% Get parameters to init server
argv = getIntMacrosFromCHeader(headerParametersFilepath,...
    {'SERVER_PORT'; 'CONSTELLATION_N_PLANES'; 'CONSTELLATION_N_PER_PLANE'});
% Setup parammeters
N = argv(2)*argv(3);
port = argv(1);
addr = '127.0.0.1';
% Setup server
% Output as vector: 0 | Output as matrix: 1
tudat = tudatMatlabServer(port,addr,N,0);
% Setup cleanup
tudatCleanUp = onCleanup(@()tudat.delete());
% Wait for tudat-app
tudat.waitForClient();
end

%% Define constellation
fprintf('@MATLAB server: Defining nominal constellation.\n');
if uploadHeaderParameters
% Define nominal constellation
argv = getIntMacrosFromCHeader(headerParametersFilepath,...
    {'CONSTELLATION_N_PLANES'; 'CONSTELLATION_N_PER_PLANE';...
     'CONSTELLATION_F'});
numberOfPlanes = argv(1);
numberOfSatellitesPerPlane = argv(2);
phasingParameter = argv(3);
argv = getFloatMacrosFromCHeader(headerParametersFilepath,...
    {'CONSTELLATION_INC_DEG';'CONSTELLATION_SMA';'SAT_Ct1'});
inclination = argv(1)*pi/180; %(rad)
walkerParameters = ... % i:T/P/F
    [inclination;... %(i) 
     numberOfPlanes*numberOfSatellitesPerPlane;... % (T) Number of satellites
     numberOfPlanes;... % (P)
     phasingParameter]; % (F)
semiMajorAxis = argv(2); % (m)
Ct1 = argv(3); % (N)
else
numberOfPlanes = 6;%72;
numberOfSatellitesPerPlane = 5;%22;
inclination = 53.0*pi/180; %(rad)
phasingParameter = 17;
walkerParameters = ... % i:T/P/F
    [inclination;... %(i) 
     numberOfPlanes*numberOfSatellitesPerPlane;... % (T) Number of satellites
     numberOfPlanes;... % (P)
     phasingParameter]; % (F)
semiMajorAxis = 6921000; % (m)
Ct1 = 0.068; % (N)
end
% Number of satellites
N = numberOfPlanes*numberOfSatellitesPerPlane;
% Init dimensions of the dynamics of each satellite
n_single = 6;
m_single = 3;
o_single_self = 4;
o_single_rel = 2;

%% Define simulation parameters
fprintf('@MATLAB server: Retrieving simulation parameters.\n');
if uploadHeaderParameters
    argv = getFloatMacrosFromCHeader(headerParametersFilepath,...
        {'EPOCH_CONTROL_UPDATE','EPOCH_START','EPOCH_END'});
    Tctrl = argv(1); % Control cycle time (s)
    ItSim = (argv(3)-argv(2))/Tctrl+1;
    if ~tudatSimulation
        if exist('matlab_check','var') && matlab_check
            Tsim = matlab_check*Tctrl;
        else
            argv = getFloatMacrosFromCHeader(headerParametersFilepath,...
                {'EPOCH_START';'EPOCH_END'});
            Tsim = argv(2)-argv(1); % Simulation time (s)
        end
        ItSim = Tsim/Tctrl+1;
        % Load initial state
        argv = getStringMacrosFromCHeader(headerParametersFilepath,...
            {'X0_FILE_PATH'});
        x0 = readmatrix(argv{1})';
    end
else
    % Time steps
    Tctrl = 10; % Control cycle time (s)
    Tsim = 10e3; % Simulation time (s)
    ItSim = Tsim/Tctrl+1;
    % Load initial state
    x0 = readmatrix('./data/x0/x0_osculating.txt')';
    %load('./data/x0/x0_osculating.mat','x0');
end
if ~tudatSimulation
    % State evolution
    x = cell(N,1);
    u = cell(N,1);
    % Set initial state
    for i = 1:N
    %parfor I = 1:N
        x{i,1} = zeros(n_single+1,ItSim);
        u{i,1} = zeros(m_single,ItSim-1);
        x{i,1}(:,1) = x0(:,i);
    end
else
    % Satate backup if tudat simulation or MATLAB crash
    x_backup = zeros(7*N,ItSim);
end

%% Initialize controller
fprintf('@MATLAB server: Initializing thrust feedback controller.\n');
% ----- Controller parameters -----
% Tracking
trackingRange = 750e3; %5000e3; %(m)
trackingmaxInNeighbourhood = 5;
LOSRange = 5014e3; %(m)
% Scheduling 
Tt = 1;
% osculating2mean 
degree_osculating2mean = 12;
% Gain synthesis parameters
MPC_T = 100; 
MPC_d = 25;
% MPC window update
MPC_K = cell(N,MPC_d);
MPC_currentWindowGain = 0;
MPC_u = zeros(m_single,N);
% Dynamics matrices
A = cell(N,1);
B = cell(N,1);
H = cell(N,N);
Q = cell(N,1);
R = cell(N,1);
o = zeros(N,1); % Number of relative tracking outputs
% Topology matrices
Di_tau_1 = cell(N,1);
Di_tau = cell(N,1);
Di_tau_minus = cell(N,1);
% Cost 
P_kl = cell(N,1);
P_kl_prev = cell(N,1);
% Be very carefull: K_tau_1{i,1} are the gains computed in i -> that will 
% be used to compute the actuation in other agents
% On the other hand, MPC_K{i,tau+1} are the gains that are used in i to
% compute the actuation at the tau-th discrete-time insant of the window
K_tau_1 = cell(N,1);
% Start a parallel pool
if decentralized
    fprintf('@MATLAB server: ');
    parpool(35);
end

%% Request loop
fprintf('@MATLAB server: Thrust feedback controller waiting for requests.\n');
% Init simulation loop for propagation in matlab
%if ~tudatSimulation, it = 0; end
it = 0;
while(1)    
    %% Init control feedback for this instant
    % Increase iteration counter
    it = it + 1;
    % Wait for actuation request
    if tudatSimulation
        [t,x_t,state] = tudat.getRequest();
        if(state),break; end
        x_backup(:,it) = x_t;
        if rem(it,MPC_d) == 0
            writematrix(x_backup,'./output/tmp_matlab_x_t.txt','Delimiter','tab');
        end
    else      
        % Check termination
        if it == ItSim, break; end
        % Prepare pseudo-input for controller feedback
        x_t = zeros((n_single+1)*N,1);
        for i = 1:N
            x_t((i-1)*(n_single+1)+1:i*(n_single+1),1) = x{i,1}(:,it);
        end
        t = (it-1)*Tctrl;
    end
    
    %% Run an iteration of the DDLQR algorithm
    DDELQR_iteration_routine;

    %% Vectorize and send constellation thrust command
    % Send actuation in TNW frame (X: velocity; Y: - along radius; Z: angular momentum)
    if tudatSimulation
        % Tudat's TNW frame is broken! Sencond and third axis are *(-1),
        % i.e., in tudat (X: velocity; Y:  along radius; Z: - angular momentum)
        u_tudat = MPC_u;
        u_tudat([2 3],:) = -u_tudat([2 3],:);
        tudat.sendResponse(u_tudat(:)); 
    else
        % Simulation 
        % From pseudo-output to structure 
        for i = 1:N
            u{i,1}(:,it) = MPC_u(:,i);%u_t((i-1)*m_single+1:i*m_single);
        end
        % Propagate orbit
        %for i = 1:N
        parfor i = 1:N
            x{i,1}(:,it+1) = propagateSatellite(x{i,1}(:,it),u{i,1}(:,it),Tctrl);
        end
    end
end
toc;

%% Save data if simulatio was performed in matlab
if ~tudatSimulation && (~exist('matlab_check','var') || ~matlab_check)
    mkdir('output');
    save('./output/output_matlab_propagation.mat','x','u','uploadHeaderParameters');
end

%% Termination controller
if tudatSimulation
% Shut down parallel pool
delete(gcp('nocreate'));
clear tudatCleanUp;
fprintf('@MATLAB server: Thrust feedback controller has been terminated.\n');
end
