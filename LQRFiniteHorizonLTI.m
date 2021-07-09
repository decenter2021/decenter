function [Kinf,Pinf] = LQRFiniteHorizonLTI(A,B,Q,R,E,opts)
%% Description
% This function computes the steady-state finite-horizon LQR regulator gain
% Input:    - A,B,Q,R
%           - E: sparsity pattern
%           - opts: optional input arguments
%               - epsl: minimum relative improvement on the objective function
%               - maxOLIt: maximum number of iterations until convergence 
%               - W: starting window length
%               - maxOLIt: maximum number of outer loop iterations
%               - findWindowLength: if enabled iterates through window 
%               length values until convergence is reached
%               - verbose: display algorithm status messages
% Output:   - Kinf: nxo steady-state gain matrix
%           - Pinf: nxn steady-state estimation error covariance matrix
% Important notes: 
%           - output gain corresponds to the control law: u(k)=-K(k)*x(k)
% WARNING: Returns Kinf = NaN and Pinf = NaN if convergence could not be reached

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

LogicalStr = {'false', 'true'};
if opts.verbose
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Running finite-horizon algorithm with:\nepsl = %g | W = %d | maxOLIt = %d | findWindowSize = %s.\n',opts.epsl,opts.W,opts.maxOLIt,LogicalStr{opts.findWindowLength+1});
end

%% Gains computation    
n = size(A,1); % Get value of n from the size of A 
% Initialise Finite Horizon with One Step gain and covariance matrices
[K,P] = LQROneStepSequenceLTI(A,B,Q,R,E,opts.W);
Z = vectorZ(vec(E)); % Compute matrix Z
Pprev = zeros(n,n); % Previous iteration
Kinf = NaN;
counterSteadyState = 0; % Counter for the number of iterations for which a steady-state solution was found
counterOuterLoopIt = 0; % Counter for the number of outer loop iterations
while  true % Outer loop
    counterOuterLoopIt = counterOuterLoopIt+1; % Increment number of outer loop iterations
    % Inner loop
    for k = opts.W:-1:1
        % Update covariance matrix after the update step
        if k > 1
            S = B'*P{k-1,1}*B+R;
        else
            S = B'*Q*B+R;
        end
        % Compute summation to obtain matrix Lambda
        Lambda = eye(n);
        for i = k:opts.W
           % Compute summation to  obtain matrix Gamma
           Gamma = eye(n);
           for j = k+1:i
              Gamma = Gamma*(A-B*K{i}); 
           end
           Lambda = Lambda + Gamma*Gamma';
        end   
        if k > 1
            C_eq = B'*P{k-1,1}*A*Lambda;
        else
            C_eq = B'*Q*A*Lambda;
        end
        % Adjust gain using efficient solver [1]
        K{k,1} = sparseEqSolver(S,Lambda,...
            C_eq,E);

    end
    % Recompute covariances 
    for i = 1:opts.W
        if i == 1
            P{i,1} = Q + K{i,1}'*R*K{i,1}+(A-B*K{i,1})'*Q*(A-B*K{i,1});
        else
            P{i,1} = Q + K{i,1}'*R*K{i,1}+(A-B*K{i,1})'*P{i-1}*(A-B*K{i,1});
        end
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
            [K,P] = LQROneStepSequenceLTI(A,B,Q,R,E,opts.W);
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
function [K,P] = LQROneStepSequenceLTI(A,B,Q,R,E,w)
%% Gain computation
n = size(A,1); % Get value of n from the size of A 
m = size(B,2); % Get value of o from the size of C 
K = cell(w,1);
P = cell(w,1);
for l = 1:w
    % Update the covariance of the update step, P_. P is the covariance 
    % matrix after the filtering step.
    if l == 1
        S = R+B'*Q*B;
    else
        S = R+B'*P{l-1,1}*B;
    end
    % Initialise gain matrix 
    K{l,1} = zeros(m,n);
    % Summation to compute the gain
    for i = 1:n
        L = zeros(n);
        L (i,i) = 1; % Generate matrix L_i
        M = zeros(m);
        for j = 1:m % Gererate matrix M_i
            if E(j,i) ~= 0
                M(j,j) = 1;
            end
        end
        % Compute the ith term of the summation 
        if l == 1
            K{l,1} = K{l,1} + (eye(m)-M+M*S*M)\M'*B'*Q*A*L';
        else
            K{l,1} = K{l,1} + (eye(m)-M+M*S*M)\M'*B'*P{l-1,1}*A*L';
        end
        
    end
  
    % Update the covariance matrix after the filtering step
    if l == 1
       P{l,1} = Q + K{l,1}'*R*K{l,1}+(A-B*K{l,1})'*Q*(A-B*K{l,1});
    else
        P{l,1} = Q + K{l,1}'*R*K{l,1}+(A-B*K{l,1})'*P{l-1,1}*(A-B*K{l,1});
    end

end
end

%[1] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497

