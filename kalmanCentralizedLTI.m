function [Kinf,Pinf] = kalmanCentralizedLTI(A,C,Q,R,opts)
%% Description
% This function computes the steady-state centralized Kalman filter gain 
% Input:    - A,C,Q,R
%           - opts: optional input arguments
%               - epsl: minimum relative improvement on the objective function
%               - maxIt: maximum number of iterations until convergence 
%               - P0: initialization estimation error covariance matrix
%               - verbose: display algorithm status messages           
% Output:   - Kinf: nxo steady-state gain matrix
%           - Pinf: nxn steady-state estimation error covariance matrix
% WARNING: Returns Kinf = NaN and Pinf = NaN if convergence could not be reached

%% Argument handling
if ~exist('opts','var') 
    opts.verbose = false; % Default is not to display algorithm status messages
elseif ~isfield(opts,'verbose')
    opts.verbose = false; % Default is not to display algorithm status messages
end
if ~isfield(opts,'maxIt')
    opts.maxIt = 1000; % Default maximum number of iterations until convergence
end
if ~isfield(opts,'epsl')
    opts.epsl = 1e-5; % Default minimum relative improvement on the objective function
end
if ~isfield(opts,'P0')
    opts.P0 = zeros(size(A)); % Default initialization estimation error covariance matrix chosen to be the null nxn matrix
end
if opts.verbose
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Computing centralized kalman filter with: epsl = %g | maxIt = %d.\n',opts.epsl,opts.maxIt);
end

%% Gain computation
n = size(A,1); % Get value of n from the size of A 
%o = size(C,1); % Get value of o from the size of C
Pprev = nan; % Previous iteration
Pinf = opts.P0;
for i = 1:opts.maxIt
    % Generate gain and update covariance according to the centralized 
    % Kalman filter update and filtering iterations. P_ is the covariance
    % matrix of the update step and P the covariance matrix after the
    % filtering step.
    P_ = A*Pinf*transpose(A)+Q;
    Kinf = P_*transpose(C)/(C*P_*transpose(C)+R);
    Pinf = Kinf*R*transpose(Kinf)+(eye(n)-Kinf*C)*P_*transpose(eye(n)-Kinf*C);

    % Check if new iteration is within the tolerance
    if  abs(trace(Pinf)-trace(Pprev))/trace(Pprev) < opts.epsl   
        if opts.verbose
            fprintf("Convergence reached with: epsl = %g | maxIt = %d.\n",opts.epsl,opts.maxIt);
            fprintf('A total of %d iterations were run.\n',i);
            fprintf('----------------------------------------------------------------------------------\n');
        end
        break; 
    elseif i == opts.maxIt % Convergence could not be reached
        fprintf("Centralized algorithm was unable to reach convergence with the specified\nparameters: epsl = %g | maxIt = %d.\n",opts.epsl,opts.maxIt);
        fprintf("Sugested actions:\n- Increase \'maxIt\', the maximum number of iterations;\n- Increase \'epsl\', the minimum relative improvement on the objective function\noptimization problem.\n");
        fprintf('----------------------------------------------------------------------------------\n');
        Kinf = NaN;
        Pinf = NaN;
    end
    Pprev = Pinf;
end

end
