%% Package: tudat-matlab-thrust-feedback
% Author: Leonardo Pedroso
%% Simulation of constellation of LEO satellites 
% Numeric dynamics simulation in tudat 
% Thrust feedback in MATLAB using package tudat-matlab-thrust-feedback
%% Init 
% Variable 'matlab_check' may be received as an input to run a small
% simulation propagated in matlab to make sure there are no errors
clearvars -except matlab_check; % Clear workspace variables
%% Add source of packages to matlab path
% Add source of tudat feedback matlab server class
addpath('/mnt/nfs/home/lpedroso/lib/src-tudat-matlab-thrust-feedback'); 
tic; % Log the execution time

%% Define execution options 
% Turn off tudat simulation for debug
tudatSimulation = true;
% Turn off parameter upload from C header
uploadHeaderParameters = true;
headerParametersFilepath = 'tudat-matlab-parameters.h';

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
numberOfPlanes = 6;
numberOfSatellitesPerPlane = 5;
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
    x0 = readmatrix('./data/x0/x0_constellation_30.txt')';
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
    % Simulation environment 
 	env.mu = 3.986004418e14; %(m^3 s^-2)
    env.RE = 6371e3; %(m)
    env.J2 = 1082.6267e-6;
    env.wE = 7.2921150e-5; %(rad/s)
    % Atmosphere cf. Table 7-4 Vallado1997
    env.h0 = 500e3; %(m)
    env.rho0 = 6.967e-13; %(Kg/m^3)
    env.H = 63.822e3; %(m)
    % Satellite
    argv = getFloatMacrosFromCHeader(headerParametersFilepath,...
    {'SAT_Cd';'SAT_Ad'});
    env.Cd = argv(1);
    env.Ad = argv(2); % (m^2)
    % Thrust
    argv = getFloatMacrosFromCHeader(headerParametersFilepath,...
    {'SAT_Ct1';'SAT_ISP';'SAT_g0'});
    env.Ct1 =  argv(1); % (N)
    env.Isp=  argv(2); % (s)
    env.g0 =  argv(3); % (ms^-2)
    env.fs = 1;
else
    % Satate backup if tudat simulation or MATLAB crash
    x_backup = zeros(7*N,ItSim);
end

%% Initialize controller
fprintf('@MATLAB server: Initializing thrust feedback controller.\n');
% ----- Controller parameters -----
% Define controller variables here
% E.g. define kp, kd, ki of a PID controller
% E.g. define window length of MPC scheme 
thrust_mag = 0.06; %(N)

%% Request loop
fprintf('@MATLAB server: Thrust feedback controller waiting for requests.\n');
% Init simulation loop for propagation in matlab
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
        if rem(it,100) == 0
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
        % Progress bar
        barWidth = 70;
        progress = it/(ItSim-1);
        fprintf("[");
        pos = round(barWidth * progress);
        for i = 0:barWidth-1
        	if i < pos, fprintf("=");
        	elseif i == pos, fprintf(">");
        	else fprintf(" ");
        	end
        end
        fprintf("] %.2f %%",progress*100);
        if progress ~= 1, fprintf("\r"); 
        else, fprintf("\n"); end
    end
    
    %% Run feedback routine
    matlab_feedback_routine;

    %% Vectorize and send constellation thrust command
    % Send actuation in TNW frame (X: velocity; Y: - along radius; Z: angular momentum)
    if tudatSimulation
        % Tudat's TNW frame is broken! Sencond and third axis are *(-1),
        % i.e., in tudat (X: velocity; Y:  along radius; Z: - angular momentum)
        u_tudat = u_t;
        u_tudat([2 3],:) = -u_tudat([2 3],:);
        tudat.sendResponse(u_tudat(:)); 
    else
        % Simulation 
        % From pseudo-output to structure 
        for i = 1:N
            u{i,1}(:,it) = u_t((i-1)*m_single+1:i*m_single);
        end
        % Propagate orbit
        for i = 1:N
        %parfor i = 1:N
            x{i,1}(:,it+1) = propagateSatellite(x{i,1}(:,it),u{i,1}(:,it),Tctrl,env);
        end
    end
end
toc;

%% Save data if simulation was performed in matlab
if ~tudatSimulation && (~exist('matlab_check','var') || ~matlab_check)
    mkdir('output');
    save('./output/output_matlab.mat','x','u','uploadHeaderParameters');
end

%% Termination controller
if tudatSimulation
% Shut down any parallel pool
delete(gcp('nocreate'));
clear tudatCleanUp;
fprintf('@MATLAB server: Thrust feedback controller has been terminated.\n');
end

