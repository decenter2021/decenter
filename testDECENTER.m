function testDECENTER()
    version = '1.0.0';
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Testing installation of DECENTER v%s.\n',version);
    currentDir = pwd;
    cd ~; % Go to the root folder to check if installation was successfull
    try
        fprintf('----------------------------------------------------------------------------------\n');
        % Check if all the functions are operational
        fprintf('kalmanCentralizedLTI: ');
        kalmanCentralizedLTI(1,1,1,1);
        fprintf('kalmanOneStepLTI: ');
        kalmanOneStepLTI(1,1,1,1,1);
        fprintf('kalmanFiniteHorizonLTI: ');
        kalmanFiniteHorizonLTI(1,1,1,1,1);
        fprintf('----------------------------------------------------------------------------------\n');
        fprintf('DECENTER v%s was successfully installed.\n',version);
    catch me
        fprintf('Installation of DECENTER v%s failed.\n',version);
        switch me.identifier
            case 'MATLAB:UndefinedFunction'
                fprintf('ERROR: DECENTER-toolbox-%s is not on the MATLAB path.\n',version);
                fprintf('Sugested action: add DECENTER-toolbox-%s to your MATLAB path on ´Set Path´under\nthe MATLAB ´Home´tab.\n',version);
            otherwise
                fprintf('ERROR: An unknown error has occured.\n');
                fprintf('Sugested action: reinstall DECENTER.\n');
        end
    end 
    fprintf('----------------------------------------------------------------------------------\n');
    cd(currentDir);
end

