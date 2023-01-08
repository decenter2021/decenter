%% Covert Tudat .dat ouput files to .mat files
%% Init workspace 
clear;
clc;
%% Get port and number of satellites
args = {'CONSTELLATION_N_PLANES', 'CONSTELLATION_N_PER_PLANE','CONSTELLATION_F',...
    'CONSTELLATION_SMA',...
    'EPOCH_START','CONSTELLATION_ECC','CONSTELLATION_AOP',...
    'EPOCH_END', 'EPOCH_SAMPLE', 'EPOCH_CONTROL_UPDATE',...
    'SAT_Cd','SAT_Ad', 'SAT_MASS', 'SAT_Ct1', 'SAT_ISP', 'SAT_g0',...
    'SAT_Cr', 'SAT_SRPA'};
argv = cell(1,length(args));
fid = fopen('simulationParameters.h');
while 1
    tline = fgetl(fid);
    for i = 1:length(args)
        len = length(args{i})+8;
        if length(tline) > len && strcmp(tline(1:len), sprintf('#define %s',args{i}))
            argv{i} = sscanf(tline,sprintf('#define %s %%d',args{i}));
        end
    end
    if ~ischar(tline), break, end % EOF 
end
fclose(fid);
params = cell2struct(argv,args,2);

%% Set parameters
NSats = argv{1}*argv{2};

%% Init variables
x = cell(NSats,1);
%% Start converting output
if NSats > 100
    % Start a parallel pool
    %parpool(24);
    %parfor i = 1:NSats
    for i = 1:NSats
        fpath = sprintf('./output/stateSat%d.dat',i-1);
        data = dlmread(fpath);
        x{i,1} = data(:,2:end)';  
        fpath = sprintf('./output/inputSat%d.dat',i-1);
        data = dlmread(fpath);
        u{i,1} = data(:,2:end)';
    end
    % Shut down parallel pool
    delete(gcp('nocreate'));
else
        for i = 1:NSats
        fpath = sprintf('./output/stateSat%d.dat',filename,i-1);
        data = dlmread(fpath);
        x{i,1} = data(:,2:end)';  
        fpath = sprintf('./output/inputSat%d.dat',filename,i-1);
        data = dlmread(fpath);
        u{i,1} = data(:,2:end)';
    end
end

%% Save cell
save("./output/output.mat",'x','u','params');
clear;
