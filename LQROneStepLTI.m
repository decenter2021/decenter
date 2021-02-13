function [Kinf,Pinf] = LQROneStepLTI(A,B,Q,R,E,opts)
%% Description
% This function computes the steady-state one-step LQR regulator gain
% Input:    - A,B,Q,R
%           - E: sparsity pattern
%           - opts: optional input arguments
%               - epsl: minimum relative improvement on the objective function
%               - maxIt: maximum number of iterations until convergence 
%               - verbose: display algorithm status messages
% Output:   - Kinf: nxo steady-state gain matrix
%           - Pinf: nxn steady-state estimation error covariance matrix
% Important notes: 
%           - output gain corresponds to the control law: u(k)=-K(k)*x(k)
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
if opts.verbose
    fprintf('----------------------------------------------------------------------------------\n');
    fprintf('Running one-step algorithm with: epsl = %g | maxIt = %d.\n',opts.epsl,opts.maxIt);
end

%% Gain computation
n = size(A,1); % Get value of n from the size of A 
m = size(B,2); % Get value of n from the size of B 
Pinf = Q; % terminal condition
Pprev = NaN; % Previous iteration
it = opts.maxIt;
while it > 0 % LQ iterations
    Kinf = zeros(m,n);
    S = R+B'*Pinf*B;
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
        Kinf = Kinf + (eye(m)-M+M*S*M)\M'*B'*Pinf*A*L';
    end
    % Update P
    Pinf = Q + Kinf'*R*Kinf+(A-B*Kinf)'*Pinf*(A-B*Kinf);
    it = it-1;
    if abs(trace(Pinf)-trace(Pprev))/trace(Pprev) < opts.epsl
        if opts.verbose
                fprintf("Convergence reached with: epsl = %g | maxIt = %d.\n",opts.epsl,opts.maxIt);
                fprintf('A total of %d iterations were run.\n',opts.maxIt-it);
                fprintf('----------------------------------------------------------------------------------\n');
        end
        break; 
    end
    Pprev = Pinf;
    if it == 0
        fprintf("One-step did not converge.\n");
        Pinf = NaN;
        Kinf = NaN;
    end
end    
end
