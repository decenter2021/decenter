function [K,P_inf] = MovingFiniteHorizonGainLTI(A,C,Q,R,E,w_ss,epsl_inf)
    % Set epsl
    epsl = epsl_inf/10;
    % Init covariance and gains over the window 
    P = cell(w_ss,1);
    K = cell(w_ss,1);
    % Null initialization
    P0 = zeros(size(A));
    P_prevIt = NaN;
    % Iterate until convergence
    k = 0;
    while true
        k = k+1;
        % If k<w_ss choose maximum possible window
        w = min(k,w_ss);
        % w < w_ss
        if k <= w_ss
            [K(1:w,1),P{k,1}] = MovingFiniteHorizonGain(A,C,Q,R,E,w,epsl,P0,1);
            if k == w_ss
                P_prevIt = P{end,1};
            end
        % w >= w_ss
        else
            % Get covariance at the end of the window
            Paux = P{1,1};
            % Shift window covariances
            for i = 1:w-1
                P{i,1} = P{i+1,1};
            end
            % Compute new moving finite horizon gains
            [K(1:w,1),P{end,1}] = MovingFiniteHorizonGain(A,C,Q,R,E,w,epsl,Paux,0,K);
            % Check convergence
            if abs(trace(P{end,1}-P_prevIt)/trace(P_prevIt)) < epsl_inf
                break; 
            end
            P_prevIt = P{end,1};
        end 
    end
    % Output projected steady-state covariance
    P_inf = P{end,1};
end

function [K,P_end] = MovingFiniteHorizonGain(A,C,Q,R,E,w,epsl,P0,flag,K)
    n = size(A,2); % Get value of n from the size of A 
    P = cell(w,1); % Initialise cell to hold all covariance matrices
    if flag % if k<W_ss
        % Initialise with One Step gain and covariance matrices
        K = cell(w,1); % Initialise cell to hold all gain matrices
        for i = 1:w
            if i == 1
                [K{i,1},P{i,1}] = OneStepGain(A,C,Q,R,E,P0);
            else
                [K{i,1},P{i,1}] = OneStepGain(A,C,Q,R,E,P{i-1,1});
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
            P{i,1} = K{i,1}*R*K{i,1}'+(eye(n)-K{i,1}*C)*P_*(eye(n)-K{i,1}*C)';
        end
    end
    
    % Outer loop
    itCount = 0;
    while true
        itCount = itCount +1;
        P_prevIt = P{end,1};
        % Inner loop
        for i = w:-1:1
            % Update covariance matrix after the update step
            if i > 1 
                P_ = A*P{i-1,1}*transpose(A)+Q;
            else
                P_ = A*P0*transpose(A)+Q;
            end
            % Compute gamma(i+1,w)
            gamma = eye(n);
            for j = i+1:w
                gamma = gamma*(eye(n)-K{i+1+w-j,1}*C)*A; 
            end
            % Compute Psi
            Psi = transpose(gamma)*gamma;
            % Adjust gain using efficient solver [2]
            K{i,1} = sparseEqSolver(Psi,C*P_*C'+R,Psi*P_*C',E);
            % Old solver commented 
            %K{i,1} = unvec(transpose(Z)/(Z*(kron(C*P_*transpose(C)+R,Psi))...
            %   *transpose(Z))*Z*vec(Psi*P_*transpose(C)),n);
        end     
       
        % Recompute covariances 
        for i = 1:w
            if i >1
                P_ = A*P{i-1,1}*A'+Q;
            else
               P_ = A*P0*A'+Q;
           end
            P{i,1} = K{i,1}*R*K{i,1}'+(eye(n)-K{i,1}*C)*P_*(eye(n)-K{i,1}*C)';
       end
       
       % Check if minumum realative improvement has been reached
       if abs(trace(P{end,1}-P_prevIt)/trace(P_prevIt)) < epsl
          break; 
       end     
    end
    % Output last covariance matrix
    P_end = P{end,1};
    
end

%% Auxiliary functions

% This function computes the One Step gain subject to a sparsity constraint 
% for the instant when it is called
% Input:    - A,C,Q,R matrices
%           - E: a matrix that defines the sparsity pattern
%           - Pprev: nxn covariance matrix of the previous instant
% Output:   - K: nxn gain matrix for the instant of the call.
%           - Pcurr: nxn covariance matrix after both the update and
%             filtering steps (this output is also used for the 
%             initialisation of the Finite Horizon algorithm)
function [K,Pcurr] = OneStepGain(A,C,Q,R,E,Pprev)
    % Gain computation
    n = size(C,2); % Get value of n from the size of C
    %o = size(C,1); % Get value of o from the size of C 
    % Update the covariance of the update step, P_. P is the covariance 
    % matrix after the filtering step.
    P_ = A*Pprev*A'+Q;
    % Compute gain with sparse matrix solver [2]
    K = sparseEqSolver(eye(n),C*P_*C'+R,P_*C',E);
    % Update the covariance matrix after the filtering step
    Pcurr = K*R*K'+...
        (eye(n)-K*C)*P_*(eye(n)-K*C)';
end

%% References 
% [2] Pedroso, Leonardo, and Pedro Batista. 2021. "Efficient Algorithm for the 
% Computation of the Solution to a Sparse Matrix Equation in Distributed Control 
% Theory" Mathematics 9, no. 13: 1497. https://doi.org/10.3390/math9131497
