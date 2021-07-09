function [Kinf,Pinf] = kalmanOneStepLTI(A,C,Q,R,E,opts)
%% Description
% This function computes the steady-state distributed one-step Kalman 
% filter gain 
% Input:    - A,C,Q,R
%           - E: a matrix that defines the sparsity pattern
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
    fprintf('Running one-step algorithm with: epsl = %g | maxIt = %d.\n',opts.epsl,opts.maxIt);
end

%% Gain computation
n = size(A,1); % Get value of n from the size of A 
Pprev = zeros(n,n); % Previous iteration
Pinf = opts.P0;
for l = 1:opts.maxIt
    % Update the covariance of the update step, P_. P is the covariance 
    % matrix after the filtering step.
    P_ = A*Pinf*transpose(A)+Q;
    % Compute gain matrix 
    Kinf = sparseEqSolver(eye(n),C*P_*transpose(C)+R,P_*transpose(C),E);  
    % Update the covariance matrix after the filtering step
    Pinf = Kinf*R*transpose(Kinf)+...
        (eye(n)-Kinf*C)*P_*transpose(eye(n)-Kinf*C);        
    % Check if new iteration is within the relative minimum improvement
    if abs(trace(Pinf-Pprev))/trace(Pprev)<opts.epsl
       if opts.verbose
           fprintf("Convergence reached with: epsl = %g | maxIt = %d.\n",opts.epsl,opts.maxIt);
           fprintf('A total of %d iterations were run.\n',l);
           fprintf('----------------------------------------------------------------------------------\n');
       end
       break; 
    elseif l == opts.maxIt % Convergence could not be reached
       fprintf("One-step algorithm was unable to reach convergence with the specified parameters:\nepsl = %g | maxIt = %d\n",opts.epsl,opts.maxIt);
       fprintf("Sugested actions:\n- Increase \'maxIt\', the maximum number of iterations;\n- Increase \'epsl\', the minimum relative improvement on the objective function\noptimization problem.\n");
       fprintf('----------------------------------------------------------------------------------\n');
       Kinf = NaN;
       Pinf = NaN;
    end
    Pprev = Pinf;
end
end

