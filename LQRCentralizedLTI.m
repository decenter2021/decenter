function [Kinf,Pinf] = LQRCentralizedLTI(A,B,Q,R,opts)
%% Description
% This function computes the steady-state centralized LQR gain 
% Input:    - A,B,Q,R
% Output:   - Kinf: nxo steady-state gain matrix
%           - Pinf: nxn steady-state matrix P
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
    fprintf('Computing centralized LQR gain with: epsl = %g | maxIt = %d.\n',opts.epsl,opts.maxIt);
end

%% Gain computation
    n = size(A,1); % Get value of n from the size of A   
    Pinf = Q; % terminal condition 
    % Generate gain and update covariance for centralized LQR
    Pprev = NaN; % Previous iteration
    it = opts.maxIt;
    while it > 0 % LQ iterations
        S = Pinf+Q;
        P_ = A'*S*A;
        Kinf = (R+B'*S*B)\B'*S*A;
        Pinf = Kinf'*R*Kinf+(eye(n)-Kinf'*B'/A')*P_*(eye(n)-Kinf'*B'/A')';
        it = it-1;
        if (trace(Pinf)-trace(Pprev))/trace(Pprev) < opts.epsl
            if opts.verbose
                fprintf("Convergence reached with: epsl = %g | maxIt = %d.\n",opts.epsl,opts.maxIt);
                fprintf('A total of %d iterations were run.\n',opts.maxIt-it);
                fprintf('----------------------------------------------------------------------------------\n');
            end
            break; 
        end
        Pprev = Pinf;
        if it == 0
            fprintf("Centralized LQR did not converge.\n");
            Pinf = NaN;
            Kinf = NaN;
        end
    end
end