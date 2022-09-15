function [Kinf,Pinf,Pseq] = MHEMovingFiniteHorizonLTI(A,C,Q,R,E,W,opts)
%% Description
% This function computes the steady state decentralized moving finite
% horizon steady-state sequence of gains according to [1]
% Input:    - A,C,Q,R
%           - E: a matrix that defines the sparsity pattern
%           - W: window length
%           - opts: optional input arguments
%               - epsl_inf: minimum relative improvement on the objective
%               function of (12) [default = 1e-5]
%               - epsl: minimum relative improvement on the objective
%               function of (16) [default = epsl_inf/10]
%               - maxIt: maximum number of iterations in Table 2. until convergence 
%               - P0: initialization estimation error covariance matrix
%               - verbose: display algorithm status messages
% Output:   - Kinf: Wx1 cell of nxo steady-state sequence of gain matrices
%           - Pinf: steady state covariance matrix
%           - Pseq: (W+1)x1 cell of nxn covariance matrices throught the window in
%           the last iteration
% Returns Kinf = NaN and Pinf = NaN if convergence could not be reached

%% Argument handling
if ~exist('opts','var') 
    opts.verbose = false; % Default is not to display algorithm status messages
elseif ~isfield(opts,'verbose')
    opts.verbose = false; % Default is not to display algorithm status messages
end
if ~isfield(opts,'maxIt')
    opts.maxIt = 100; % Default maximum number of iterations until convergence
end
if ~isfield(opts,'epsl_inf')
    opts.epsl_inf = 1e-4; % Default minimum relative improvement on the objective function
end
if ~isfield(opts,'epsl')
    opts.epsl = opts.epsl_inf/10; % Default minimum relative improvement on the objective function
end
if ~isfield(opts,'P0')
    opts.P0 = zeros(size(A)); % Default initialization estimation error covariance matrix chosen to be the null nxn matrix
end
if opts.verbose
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Running moving finite horizon algorithm with:\nepsl_inf = %g | W = %d | maxIt = %d .\n',opts.epsl_inf,W,opts.maxIt);
end

%% Gains computation

% Initialize gain and covarinace over the moving window
P = cell(W,1);
Pseq = cell(W+1,1);
Kinf = cell(W,1);

% Table 2. in [1]
k = 0;
P_prevIt = NaN;
while true
    k = k+1;
    % Select maximum window length while k<= W
    w = min(k,W);
    if k <= W
        % Compute new sequence
        [Kinf(1:w,1),P{k,1}] = MovingFiniteHorizonGain(A,C,Q,R,E,w,opts.epsl,opts.P0,1);
        if k == W
            P_prevIt = P{end,1};
        end
    else
        Pseq{1,1} = P{1,1};
        for i = 1:w-1
            P{i,1} = P{i+1,1};
        end
        % Compute new sequence
        [Kinf(1:w,1),P{end,1},Pseq(2:W+1,1)] = MovingFiniteHorizonGain(A,C,Q,R,E,w,opts.epsl,Pseq{1,1},0,Kinf);    
        % Check convergence
        if abs(trace(P{end,1}-P_prevIt)/trace(P_prevIt)) < opts.epsl_inf
            if opts.verbose
                fprintf(sprintf("Convergence reached with: epsl_inf = %g | W = %d | maxOLIt = %d\n",opts.epsl_inf,W,opts.maxIt));
                fprintf('A total of %d outer loop iterations were run.\n',k);
                fprintf('----------------------------------------------------------------------------------\n');
            end
            break; 
        elseif k == opts.maxIt
            fprintf("Moving Finite-horizon algorithm was unable to reach convergence with the specified\nparameters: epsl_inf = %g | W = %d | maxIt = %d\n",opts.epsl_inf,W,opts.maxIt);
            fprintf("A total of %d outer loop iterations were run.\n",k);
            fprintf("Sugested actions:\n- Manually tune \'epsl_inf' and \'maxIt\' (in this order);\n- Increase \'epsl_inf\', the minimum relative improvement on the objective function\noptimization problem.\n- Increase \'maxIt\', the maximum number of iterations.\n");
            fprintf('----------------------------------------------------------------------------------\n');
            Kinf = NaN;
            Pinf = NaN;
            break;
        end
        P_prevIt = P{end,1};
    end 
end
% Output steady-state covariance
Pinf = P{end,1};

end

%% Compute a window of moving finite horizon gains 
% Algorithm in Table 1 [1]
function [K,P_end,P] = MovingFiniteHorizonGain(A,C,Q,R,E,w,epsl,P0,flag,K)
    n = size(A,2); % Get value of n from the size of A 
    P = cell(w,1); % Initialise cell to hold all covariance matrices
    if flag % if k < W_ss
        % Initialise with One-Step gain and covariance matrices
        K = cell(w,1); % Initialise cell to hold all gain matrices
        for i = 1:w
            if i == 1
                opts.P0 = P0;
                opts.epsl = epsl;
                [K{i,1},P{i,1}] = kalmanOneStepLTI(A,C,Q,R,E,opts);
            else
                opts.P0 = P{i-1,1};
                opts.epsl = epsl;
                [K{i,1},P{i,1}] = kalmanOneStepLTI(A,C,Q,R,E,opts);
            end   
        end
    else
        % Update covariance matrices
        for i = 1:w
           if i >1
               P_ = A*P{i-1,1}*A'+Q;
           else
               P_ = A*P0*A'+Q;
           end
          P{i,1} = K{i,1}*R*K{i,1}'+(eye(n)-K{i,1}*C)*P_*((eye(n)-K{i,1}*C)');
       end
    end
    
    % Outer loop iterations counter
    itCount = 0;
    while true
        % Increse iteration counter
        itCount = itCount +1;
        P_prevIt = P{end,1};
        
        % Inner loop
        for i = w:-1:1
            % Update covariance matrix after the update step
            if i > 1 
                P_ = A*P{i-1,1}*A'+Q;
            else
                P_ = A*P0*A'+Q;
            end
            % Compute gamma(i+1,w)
            gamma = eye(n);
            for j = i+1:w
                gamma = gamma*(eye(n)-K{i+1+w-j,1}*C)*A; 
            end
            % Compute Psi
            Psi = transpose(gamma)*gamma;

            % Adjust gain using efficient solver [2]
            K{i,1} = sparseEqSolver(Psi,C*P_*C'+R,...
                Psi*P_*C',E);
        end     
        % Recompute covariances 
        for i = 1:w
            if i >1
                P_ = A*P{i-1,1}*A'+Q;
            else
                P_ = A*P0*A'+Q;
            end
            P{i,1} = K{i,1}*R*transpose(K{i,1})+(eye(n)-K{i,1}*C)*P_*transpose(eye(n)-K{i,1}*C);
        end
      
        % Check if minumum realative improvement has been reached
        if abs(trace(P{end,1}-P_prevIt)/trace(P_prevIt)) < epsl
            break; 
        end     
    end
    % Output last covariance matrix
    P_end = P{end,1};  
end

%% References

% [1] Pedroso, L., and Batista, P., 2022. "Decentralized discrete-time 
% moving horizon estimation for large-scale networks of interconnected 
% unconstrained linear systems" [not published yet]

% [2] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497

