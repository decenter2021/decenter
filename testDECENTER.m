function testDECENTER()
    version = '1.2.0';
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Testing installation of DECENTER v%s.\n',version);
    currentDir = pwd;
    cd ~; % Go to the root folder to check if installation was successfull
    try
        fprintf('----------------------------------------------------------------------------------\n');
        % Check if all the functions are operational
        fprintf('kalmanCentralizedLTI: ');
        kalmanCentralizedLTI(1,1,1,1);
        fprintf('\t\tOK\n');

        fprintf('kalmanOneStepLTI: ');
        kalmanOneStepLTI(1,1,1,1,1);
        fprintf('\t\tOK\n');

        fprintf('kalmanFiniteHorizonLTI: ');
        opts = struct('W',10);
        kalmanFiniteHorizonLTI(1,1,1,1,1,opts);
        fprintf('\tOK\n');

        fprintf('LQRCentralizedLTI: ');
        LQRCentralizedLTI(1,1,1,1);
        fprintf('\t\tOK\n');

        fprintf('LQROneStepLTI: ');
        LQROneStepLTI(1,1,1,1,1);
        fprintf('\t\t\tOK\n');

        fprintf('LQRFiniteHorizonLTI: ');
        opts = struct('W',10);
        LQRFiniteHorizonLTI(1,1,1,1,1,opts);
        fprintf('\t\tOK\n');

        fprintf('----------------------------------------------------------------------------------\n');
        fprintf('DECENTER v%s was successfully installed.\n',version);
    catch me
        fprintf('Failed\n');
        fprintf('Installation of DECENTER v%s failed.\n',version);
        switch me.identifier
            case 'MATLAB:UndefinedFunction'
                fprintf('ERROR: decenter-%s is not on the MATLAB path.\n',version);
                fprintf('Sugested action: add decenter-%s to your MATLAB path on Set Path under the\nMATLAB Home tab.\n',version);
            otherwise
                fprintf('ERROR: An unknown error has occured.\n');
                fprintf('Sugested action: reinstall DECENTER.\n');
        end
    end 
    fprintf('----------------------------------------------------------------------------------\n');
    cd(currentDir);
end

