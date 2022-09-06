function [Kinf,Pinf] = kalmanFiniteHorizonLTI(A,C,Q,R,E,opts)
%% Description
% This function computes the steady state distributed finite-horizon Kalman
% filter gain according to [1]
% Input:    - A,C,Q,R
%           - E: a matrix that defines the sparsity pattern
%           - opts: optional input arguments
%               - epsl: minimum relative improvement on the objective function
%               - findWindowLength: ajust window length
%               - W: if findWindowLength = 'true' 
%                       starting value of the  window length
%                    else 
%                       single value of the finite winodw length to test
%               - maxOLIt: maximum number of outer loop iterations until convergence 
%               - P0: initialization estimation error covariance matrix
%               - verbose: display algorithm status messages
% Output:   - Kinf: nxn steady state gain matrix
%           - Pinf: steady state covariance matrix
% Returns Kinf = NaN and Pinf = NaN if convergence could not be reached

%% Argument handling verbose
if ~exist('opts','var') 
    opts.verbose = false; % Default is not to display algorithm status messages
elseif ~isfield(opts,'verbose')
    opts.verbose = false; % Default is not to display algorithm status messages
end
if ~isfield(opts,'findWindowLength')
    opts.findWindowLength = false; % Default action is to find a suitable window length
end
if ~isfield(opts,'maxOLIt')
    opts.maxOLIt = 100; % Default maximum number of iterations until convergence
end
if ~isfield(opts,'W')
    opts.W = round(2/min(abs(eig(A)))); % Default starting window size is two times the slowest natural frequency among the poles of the open-loop system 
end
if ~isfield(opts,'epsl')
    opts.epsl = 1e-5; % Default minimum relative improvement on the objective function
end
if ~isfield(opts,'P0')
    opts.P0 = zeros(size(A)); % Default initialization estimation error covariance matrix chosen to be the null nxn matrix
end

LogicalStr = {'false', 'true'};
if opts.verbose
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Running finite-horizon algorithm with:\nepsl = %g | W = %d | maxOLIt = %d | findWindowSize = %s.\n',opts.epsl,opts.W,opts.maxOLIt,LogicalStr{opts.findWindowLength+1});
end

%% Gains computation    
n = size(A,1); % Get value of n from the size of A 

% Initialise Finite Horizon with One Step gain and covariance matrices
[K,P] = OneStepSequenceLTI(A,C,Q,R,E,opts.W,opts.P0,opts.verbose);

Pprev = zeros(n,n); % Previous iteration
Kinf = NaN;
P0 = opts.P0;
counterSteadyState = 0; % Counter for the number of iterations for which a steady-state solution was found
counterOuterLoopIt = 0; % Counter for the number of outer loop iterations
while  true % Outer loop
    counterOuterLoopIt = counterOuterLoopIt+1; % Increment number of outer loop iterations
    if opts.verbose
        fprintf("Outer-loop iteration %d.\n",counterOuterLoopIt);
    end
    % Inner loop
    for i = opts.W:-1:1
        if opts.verbose
            fprintf("\tInner-loop iteration %d.\n",opts.W-i+1);
        end
        % Update covariance matrix after the update step
        if i > 1 
            P_ = A*P{i-1,1}*A'+Q;
        else
            P_ = A*P0*A'+Q;
        end
        % Compute summation to obtain matrix Lambda
        Lambda = eye(n);
        for m = i+1:opts.W
           % Compute summation to  obtain matrix Gamma
           Gamma = eye(n);
           for j = i+1:m
              Gamma = Gamma*(eye(n)-K{i+1+m-j,1}*C)*A; 
           end
           Lambda = Lambda + transpose(Gamma)*Gamma;
        end      
        % Adjust gain using efficient solver [2]
        K{i,1} = sparseEqSolver(Lambda,C*P_*transpose(C)+R,...
            Lambda*P_*transpose(C),E);
        % Old solver commented out
        % K{i,1} = unvec(transpose(Z)/(Z*(kron(C*P_*transpose(C)+R,Lambda))...
           %*transpose(Z))*Z*vec(Lambda*P_*transpose(C)),n);
    end
    % Recompute covariances 
    for i = 1:opts.W
        if i >1
            P_ = A*P{i-1,1}*A'+Q;
        else
            P_ = A*P0*A' + Q;
        end
        P{i,1} = K{i,1}*R*transpose(K{i,1})+...
          (eye(n)-K{i,1}*C)*P_*transpose(eye(n)-K{i,1}*C);
    end
    
    % Find steady-state estimation error covariance for each outer loop
    % iteration. For a precision epsl between successive outer loop
    % iterations, a precision of at most epsl/10 is required for the
    % convergence of the outer loop iteration.
    Pinf = NaN;
    rel_conv = zeros(opts.W,1);
    rel_conv(1) = nan;
    for l = 2:opts.W
        rel_conv(l) = abs(trace(P{l,1})-trace(P{l-1,1}))/trace(P{l-1,1});
    end
    [min_rel_conv, idx] = min(rel_conv);
    if min_rel_conv < opts.epsl/10
        % Increse number of times a steady-state solution was found within the
        % minimum relative improvement
        counterSteadyState = counterSteadyState+1;
        Pinf = P{idx,1};
        Kinf = K{idx,1}; 
    end

    if opts.verbose
        fprintf("Outer-loop iteration %d finished.\n",counterOuterLoopIt);
        fprintf("Window convergence within %g.\n",min_rel_conv);
        fprintf("Trace: %g\n", trace(P{end,1}));
        fprintf("LTI gain convergence within %g.\n",abs(trace(Pinf)-trace(Pprev))/trace(Pprev));
    end

    % Check if this iteration is within the tolerance
    if abs(trace(Pinf)-trace(Pprev))/trace(Pprev) < opts.epsl
        if opts.verbose
            fprintf(sprintf("Convergence reached with: epsl = %g | W = %d | maxOLIt = %d\n",opts.epsl,opts.W,opts.maxOLIt));
            fprintf('A total of %d outer loop iterations were run, out of which %.01f%% converged within\nthe specified minimum improvement.\n',counterOuterLoopIt,100*counterSteadyState/counterOuterLoopIt);
            fprintf('----------------------------------------------------------------------------------\n');
        end
        break;
    elseif counterOuterLoopIt == opts.maxOLIt  % Convergence could not be reached
        % maximum window length is 100 times the slowest natural frequency among the poles of the open-loop system 
        if ~opts.findWindowLength || opts.W > round(100/min(abs(eig(A)))) 
            fprintf("Finite-horizon algorithm was unable to reach convergence with the specified\nparameters: epsl = %g | W = %d | maxOLIt = %d\n",opts.epsl,opts.W,opts.maxOLIt);
            fprintf("A total of %d outer loop iterations were run, out of which %.01f%% converged within\nthe specified minimum improvement.\n",counterOuterLoopIt,100*counterSteadyState/counterOuterLoopIt);
            fprintf("Sugested actions:\n- Manually tune \'W\', \'epsl' and \'maxOLIt\' (in this order);\n- Increase \'W\', the finite window length;\n- Increase \'epsl\', the minimum relative improvement on the objective function\noptimization problem.\n- Increase \'maxOLIt\', the maximum number of outer loop iterations.\n");
            fprintf('----------------------------------------------------------------------------------\n');
            Kinf = NaN;
            Pinf = NaN;
            break;
        else
            counterOuterLoopIt = 0; % Restart outer loop iteration count
            opts.W = round(1.5*opts.W); % Try new window length 1.5 times bigger than the previous 
            if opts.verbose 
                fprintf('Trying new window length W = %d\n',opts.W);
            end
            % Reinitialize finite-horizon with one-step gain and covariance matrices
            [K,P] = OneStepSequenceLTI(A,C,Q,R,E,opts.W,opts.P0,opts.verbose);
        end
    end
    Pprev = Pinf;
end  

end

%% Auxiliary functions

% This function computes the sequence of Kalman filter gains for a window
% using the one-step method
% Input:    - A,C,Q,R
%           - E: a matrix that defines the sparsity pattern
%           - w: window size
% Output:   - Kinf: wx1 cell of gain matrices 
%           - Pinf: wx1 cell of estimation erro covariance matrices
function [K,P] = OneStepSequenceLTI(A,C,Q,R,E,w,P0,verbose)
%% Gain computation
n = size(A,1); % Get value of n from the size of A 
% o = size(C,1); % Get value of o from the size of C 
K = cell(w,1);
P = cell(w,1);
if verbose
    fprintf("Outer-loop initialization.\n");
end
for l = 1:w
    if verbose
        fprintf("\tOuter-loop initialization iteration: %d.\n",l);
    end
    % Update the covariance of the update step, P_. P is the covariance 
    % matrix after the filtering step.
    if l == 1
        P_ = A*P0*transpose(A)+Q;
    else
        P_ = A*P{l-1,1}*transpose(A)+Q;
    end
    % Compute gain
    K{l,1} = sparseEqSolver(eye(n),C*P_*C'+R,P_*C',E);
    % Update the covariance matrix after the filtering step
    P{l,1} = K{l,1}*R*transpose(K{l,1})+...
        (eye(n)-K{l,1}*C)*P_*transpose(eye(n)-K{l,1}*C);
end
end

%% References
% [1] Viegas, D., Batista, P., Oliveira, P. and Silvestre, C., 2018. Discrete-time 
% distributed Kalman filter design for formations of autonomous vehicles. 
% Control Engineering Practice, 75, pp.55-68.

% [2] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497

