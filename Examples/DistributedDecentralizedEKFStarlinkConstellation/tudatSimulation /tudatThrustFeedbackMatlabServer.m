%% Matlab server for thrust feedback
%% Init workspace 
clear;
clc;

%% Get port and number of satellites
args = {'SERVER_PORT'; 'CONSTELLATION_N_PLANES'; 'CONSTELLATION_N_PER_PLANE'};
argv = zeros(length(args),1);
fid = fopen('simulationParameters.h');
while 1
    tline = fgetl(fid);
    for i = 1:length(args)
        len = length(args{i})+8;
        if length(tline) > len && strcmp(tline(1:len), sprintf('#define %s',args{i}))
            argv(i) = sscanf(tline,sprintf('#define %s %%d',args{i}));
        end
    end
    if ~ischar(tline), break, end % EOF 
end
fclose(fid);

%% Initialize server
% Setup parammeters
NSats = argv(2)*argv(3);
port = argv(1);
addr = '127.0.0.1';
% Setup server
% Output as vector:0 | Output as matrix: 1
tudat = tudatMatlabServer(port,addr,NSats,1);
tudatCleanUp = onCleanup(@()tudat.delete());

% Wait for Tudat App
tudat.waitForClient();

%% Initialization controller
fprintf('@MATLAB server: Initializing thrust feedback controller.\n');
% Start a parallel pool
% parpool(24);

%% Request loop
fprintf('@MATLAB server: Thrust feedback controller waiting for requests.\n');

% Constants
Ct1 = 0.068; % (N)
Isp = 1640; % (s)
g0 = 9.81; % (ms^-2)

while(1)
    % Get state
    [t,x,state] = tudat.getRequest();
    if(state),break; end
    
    % Compute actuation  
    u = zeros(3,NSats);
      
    % Send actuation in TNW frame (X: velocity; Y: - along radius; Z: angular momentum)
    tudat.sendResponse(u(:)); 
end

%% Termination controller
% Shut down parallel pool
delete(gcp('nocreate'));
clear tudatCleanUp;
fprintf('@MATLAB server: Thrust feedback controller has been terminated.\n');

