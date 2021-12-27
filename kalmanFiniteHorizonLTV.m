function [K,P] = kalmanFiniteHorizonLTV(system,E,T,P0,opts)
%% Description
% This function computes the finite-horizon kalman filter gain matrices 
% subject to a sparsity constraint for all the instants of a window
% {k,...,k+T-1}
% Input:    - system: Tx4 cell whose rows contain matrices A,C,Q and R
%             for the whole window, i.e.,
%               - system{i,1} = A(k+i-1), i = 1,...,T
%               - system{i,2} = B(k+i-1), i = 1,...,T
%               - system{i,3} = Q(k+i-1), i = 1,...,T
%               - system{i,4} = R(k+i-1), i = 1,...,T
%           - E: a matrix that defines the sparsity pattern
%           - T: Finite Horizon time window
%           - P0 : nxn initial predicted covariance matrix
%           - opts: optional input arguments
%               - epsl: minimum relative improvement on the objective function
%               - maxOLIt: maximum number of outer loop iterations until convergence 
%               - verbose: display algorithm status messages
% Output:   - K: Tx1 cell of gain matrices for all the iterations, i.e.,
%               K(k+i-1), i = 1,...,T
%           - P: Tx1 cell of covariance matrices for all the iterations,
%           i.e., P(k+i-1), i = 1,...,T
%% Argument handling
if ~exist('opts','var') 
    opts.verbose = false; % Default is not to display algorithm status messages
elseif ~isfield(opts,'verbose')
    opts.verbose = false; % Default is not to display algorithm status messages
end
if ~isfield(opts,'maxOLIt')
    opts.maxOLIt = 100; % Default maximum number of iterations until convergence
end
if ~isfield(opts,'epsl')
    opts.epsl = 1e-5; % Default minimum relative improvement on the objective function
end
if opts.verbose
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Running finite-horizon algorithm with:\nepsl = %g | maxOLIt = %d.\n',opts.epsl,opts.maxOLIt);
end
%% Gain Computation
n = size(system{1,2},2); % Get value of n from the size of A 
K = cell(T,1); % Initialise cell to hold all gain matrices
P = cell(T,1); % Initialise cell to hold all covariance matrices
PprevIt = cell(T,1); % Initialise cell to hold all covariance matrices of the previous outer loop iteration
% Initialise Finite Horizon with one step gain and covariance matrices
for i = 1:T
    if i == 1
        [K{i,1},Ppred,P{i,1}] = kalmanOneStepLTV(system(i,:),E,P0);
    else
        [K{i,1},Ppred,P{i,1}] = kalmanOneStepLTV(system(i,:),E,Ppred);
    end
end
% Z = vectorZ(vec(E)); % Compute matrix Z
% Outer loop
for k = 1:opts.maxOLIt
    % Inner loop
    for i = T:-1:1
        % Update covariance matrix after the predict step
        if i > 1 
            P_ = system{i-1,1}*P{i-1,1}*transpose(system{i-1,1})+system{i-1,3};
        else
            P_ = P0;
        end
        % Compute summation to obtain matrix Lambda
        Lambda = eye(n);
        for m = i+1:T
            % Compute summation to  obtain matrix Gamma
            Gamma = eye(n);
            for j = i+1:m
                Gamma = Gamma*(eye(n)-K{i+1+m-j,1}*system{i+1+m-j,2})...
                 *system{i+1+m-j-1,1}; 
            end
            Lambda = Lambda + transpose(Gamma)*Gamma;
        end     
        % Adjust gain using efficient solver [1]
        K{i,1} = sparseEqSolver(Lambda,system{i,2}*P_*transpose(system{i,2})+system{i,4},Lambda*P_*transpose(system{i,2}),E);
        % Old solver commented 
        % K{i,1} = unvec(transpose(Z)/(Z*(kron(system{i,2}*P_*transpose(system{i,2})+system{i,4},Lambda))...
        %    *transpose(Z))*Z*vec(Lambda*P_*transpose(system{i,2})),n);
    end
    % Recompute covariances 
    for i = 1:T
        if i >1
            P_ = system{i-1,1}*P{i-1,1}*transpose(system{i-1,1})+system{i-1,3};
        else
            P_ = P0;
        end
        P{i,1} = K{i,1}*system{i,4}*transpose(K{i,1})+...
            (eye(n)-K{i,1}*system{i,2})*P_*transpose(eye(n)-K{i,1}*system{i,2});
    end 
    % Check convergence
    if k ~= 1
        relDif = 0;
        for l = 1:T
            relDif = max(relDif,abs(trace(P{l,1})-trace(PprevIt{l,1}))/trace(PprevIt{l,1}));
        end
        if relDif < opts.epsl
            if opts.verbose
                fprintf(sprintf("Convergence reached with: epsl = %g\n",opts.epsl));
                fprintf('A total of %d outer loop iterations were run.\n',k);
                fprintf('----------------------------------------------------------------------------------\n');
            end
            break;
        elseif k == opts.maxOLIt
            if opts.verbose
                fprintf("Finite-horizon algorithm was unable to reach convergence with the specified\nparameters: epsl = %g | T = %d | maxOLIt = %d\n",opts.epsl,T,opts.maxOLIt);
                fprintf("Sugested actions:\n- Manually tune \'epsl' and \'maxOLIt\' (in this order);\n- Increase \'epsl\', the minimum relative improvement on the objective function\noptimization problem.\n- Increase \'maxOLIt\', the maximum number of outer loop iterations.\n");
                fprintf('----------------------------------------------------------------------------------\n');
            end
        end
    end
    PprevIt = P;
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

%[1] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497
