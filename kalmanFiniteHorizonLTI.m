function [Kinf,Pinf] = kalmanFiniteHorizonLTI(A,C,Q,R,E,opts)
%% Description
% This function computes the steady state distributed finite-horizon Kalman
% filter gain 
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
    opts.findWindowLength = true; % Default action is to find a suitable window length
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
[K,P] = OneStepSequenceLTI(A,C,Q,R,E,opts.W,opts.P0);
Z = vectorZ(vec(E)); % Compute matrix Z
Pprev = zeros(n,n); % Previous iteration
Kinf = NaN;
counterSteadyState = 0; % Counter for the number of iterations for which a steady-state solution was found
counterOuterLoopIt = 0; % Counter for the number of outer loop iterations
while  true % Outer loop
    counterOuterLoopIt = counterOuterLoopIt+1; % Increment number of outer loop iterations
    % Inner loop
    for i = opts.W:-1:1
        % Update covariance matrix after the update step
        if i > 1 
            P_ = A*P{i-1,1}*transpose(A)+Q;
        else
            P_ = Q;
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
        % Adjust gain using efficient solver
        K{i,1} = sparseEqSolver(Lambda,C*P_*transpose(C)+R,...
            Lambda*P_*transpose(C),E);
        % Old solver commented out
        % K{i,1} = unvec(transpose(Z)/(Z*(kron(C*P_*transpose(C)+R,Lambda))...
           %*transpose(Z))*Z*vec(Lambda*P_*transpose(C)),n);
    end
    % Recompute covariances 
    for i = 1:opts.W
        if i >1
            P_ = A*P{i-1,1}*transpose(A)+Q;
        else
            P_ = Q;
        end
        P{i,1} = K{i,1}*R*transpose(K{i,1})+...
          (eye(n)-K{i,1}*C)*P_*transpose(eye(n)-K{i,1}*C);
    end
    % Find steady-state estimation error covariance for each outer loop
    % iteration. For a precision epsl between successive outer loop
    % iterations, a precision of at most epsl/10 is required for the
    % convergence of the outer loop iteration.
    Pinf = NaN;
    for l = 2:opts.W
        if abs(trace(P{l,1})-trace(P{l-1,1}))/trace(P{l-1,1}) < opts.epsl/10
           Kinf = K{l,1};
           Pinf = P{l,1};
           % Increse number of times a steady-state solution was found
           counterSteadyState = counterSteadyState+1; 
           break; 
        end
    end
    
    % Check if this iteration is within the tolerance
    if abs(trace(Pinf)-trace(Pprev))/trace(Pprev) < opts.epsl
        fprintf(sprintf("Convergence reached with: epsl = %g | W = %d | maxOLIt = %d\n",opts.epsl,opts.W,opts.maxOLIt));
        if opts.verbose
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
            [K,P] = OneStepSequenceLTI(A,C,Q,R,E,opts.W,opts.P0);
        end
    end
    Pprev = Pinf;
end  
end

%% Auxiliary functions
% Function that computes the vectorisation of a matrix
% Input:    - in: matrix to be vectorised
% Output:   - out: vec(in) 
function out = vec(in)
    out = zeros(size(in,2)*size(in,1),1);
    for j = 1:size(in,2)
       for i = 1:size(in,1)
           out((j-1)*size(in,1)+i) = in(i,j);
       end
    end
end

% Function that returns a matrix given its vectorisation and number of rows
% Input:    - in: vectorisation of a matrix
%           - n: number of rows of the matrix whose vectorisation is
%             input variable in
% Output:   - out: matrix with n rows whose vectorisation is input variable
%             in 
function out = unvec(in,n)
    out = zeros(n,size(in,1)/n);
    for j = 1:size(in,1)
       out(rem(j-1,n)+1, round(floor((j-1)/n))+1) = in(j); 
    end
end

% Function which computes matrix Z such that Z*vec(K) contains the non-zero
% elements of K according to the desired sparsity pattern
% Input:    - vecE: the vectorisation of the matrix that defines the
%             sparsity patern
% Output:   - Z
function Z = vectorZ(vecE)
    vecE = vecE~=0; % Normalise the sparsity pattern to a logical array
    Z = zeros(sum(vecE),size(vecE,1)); % Initialise matrix Z
    nZeros = 0;
    for i = 1:size(vecE,1)
       if vecE(i) ~= 0
           Z(i-nZeros,i) = 1; 
       else
           nZeros = nZeros+1;
       end     
    end
end


% This function computes the sequence of Kalman filter gains for a window
% using the one-step method
% Input:    - A,C,Q,R
%           - E: a matrix that defines the sparsity pattern
%           - w: window size
% Output:   - Kinf: wx1 cell of gain matrices 
%           - Pinf: wx1 cell of estimation erro covariance matrices
function [K,P] = OneStepSequenceLTI(A,C,Q,R,E,w,P0)
%% Gain computation
n = size(A,1); % Get value of n from the size of A 
o = size(C,1); % Get value of o from the size of C 
K = cell(w,1);
P = cell(w,1);
for l = 1:w
    % Update the covariance of the update step, P_. P is the covariance 
    % matrix after the filtering step.
    if l == 1
        P_ = A*P0*transpose(A)+Q;
    else
        P_ = A*P{l-1,1}*transpose(A)+Q;
    end
    % Initialise gain matrix 
    K{l,1} = zeros(n,o);
    % Summation to compute the gain
    for i = 1:n
        L = zeros(n);
        L (i,i) = 1; % Generate matrix L_i
        M = zeros(o);
        for k = 1:o % Gererate matrix M_i
           if E(i,k) ~= 0
               M(k,k) = 1;
           end
        end
        % Compute the ith term of the summation 
        K{l,1} = K{l,1} + (L*P_*transpose(C)*M)/...
            (eye(o)-M+M*(C*P_*transpose(C)+R)*M);
    end
    % Update the covariance matrix after the filtering step
    P{l,1} = K{l,1}*R*transpose(K{l,1})+...
        (eye(n)-K{l,1}*C)*P_*transpose(eye(n)-K{l,1}*C);
end
end
