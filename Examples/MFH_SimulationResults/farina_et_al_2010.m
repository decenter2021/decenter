% Farina et al. (2010) PMHE1 algorithm
% Compute estimate at t
% Inputs:
% A: global matrix A
% C: global matrix C
% Q: global matrix Q
% R: global matrix R
% B: global matrix B
% x_hat_t_1: hat{x}_{k/t-1}, k = t-N,...,t-1
% y: y^i_k, k = t-N+1,...,t
% u: u_k, k = t-N,...,t-1
% pi_t_N_1_t_2: pi^[i]_{t-N-1/t-2}
% Outputs:
% x: hat{x}^i_{k/t}, k = t-N,...,t
% pi_t_N_t_1: pi^[i]_{t-N/t-1}
function [x_hat_i_t,pi_t_N_t_1] = farina_et_al_2010(A,C,Q,R,B,x_hat_t_1,y,u,pi_i_t_N_1_t_2) 
    % Retrieve dimensions
    M = size(pi_i_t_N_1_t_2,1);
    n_i = size(A,1)/M;
    o_i = size(C,1)/M;
    N = size(x_hat_t_1,2);
    m_i = size(B,2)/M;
    % Init variables
    x_hat_i_t = zeros(n_i*M,N+1);
    x_hat_i_t_aux = cell(M,1);
    pi_t_N_t_1 = cell(M,1);
    % Compute predicted state
    x_hat_t_1_t = A*x_hat_t_1(:,end)+B*u(:,end);
    % Compute local dynamics matrices
    A_i = cell(M,1);
    C_i = cell(M,1);
    Q_i = cell(M,1);
    R_i = cell(M,1);
    for i = 1:M
        A_i{i,1} = A((i-1)*n_i+1:i*n_i,(i-1)*n_i+1:i*n_i);
        C_i{i,1} = C((i-1)*o_i+1:i*o_i,(i-1)*n_i+1:i*n_i);
        Q_i{i,1} = Q((i-1)*n_i+1:i*n_i,(i-1)*n_i+1:i*n_i);
        R_i{i,1} = R((i-1)*o_i+1:i*o_i,(i-1)*o_i+1:i*o_i);
    end
    % Distributed PMHE1
    %for i = 1:M  
    parfor i = 1:M     
        % Compute y^i_k, k = t-N+1,...,t
        y_i = y((i-1)*o_i+1:i*o_i,:);
        % Compute hat{x}^i_{t-N/t-1}
        x_hat_i_t_N_t_1 = x_hat_t_1((i-1)*n_i+1:i*n_i,1);       
        % Compute tilde{A}^[i]*hat{x}_{k/t-1}, k = t-N,...,t-1
        u_i_x = zeros(n_i,N);
        for k = 1:N
            u_i_x(:,k) = A((i-1)*n_i+1:i*n_i,:)*x_hat_t_1(:,k)-A_i{i}*x_hat_t_1((i-1)*n_i+1:i*n_i,k) + ...
                B((i-1)*n_i+1:i*n_i,(i-1)*m_i+1:i*m_i)*u((i-1)*m_i+1:i*m_i,k);
        end
        % Compute tilde{C}^[i]*hat{x}_{k/t-1}, k = t-N+1,...,t
        u_i_y = zeros(o_i,N);
        for k = 1:N-1
            u_i_y(:,k) = C((i-1)*o_i+1:i*o_i,:)*x_hat_t_1(:,k+1)-C_i{i}*x_hat_t_1((i-1)*n_i+1:i*n_i,k+1);
        end
        u_i_y(:,N) = C((i-1)*o_i+1:i*o_i,:)*x_hat_t_1_t-C_i{i}*x_hat_t_1_t((i-1)*n_i+1:i*n_i,1);
        % Compute pi_i_t_N_t_1
        pi_t_N_t_1{i,1} = propagate_pi_i(pi_i_t_N_1_t_2{i,1},A_i{i,1},C_i{i,1},Q_i{i,1},R_i{i,1});
        % Compute hat{x}^i_{k/t}, k = t-N,...,t
        x_hat_i_t_aux{i,1} = MHE_i_problem(y_i,u_i_y,u_i_x,A_i{i,1},C_i{i,1},Q_i{i,1},R_i{i,1},x_hat_i_t_N_t_1,pi_t_N_t_1{i,1});
    end
    % Concatenate global state
    for i = 1:M
         x_hat_i_t((i-1)*n_i+1:i*n_i,:) = x_hat_i_t_aux{i,1};
    end
end


% Farina et al. PMHE1 problem 
% Computation of pi^[i]_{t-N/t-1} according to (16b) of [1]
% Inputs 
% pi_t_N_1: pi^[i]_{t-N-1/t-2}
% A_i:      matrix A^[i]
% C_i:      matrix C^[i]
% Q_i:      process noise covariance estimated in the last time instants,
%           i.e., Q^i_{k/t-2}, k = t-N-1,...,t-2
%           USING APPROACH B: input is a single matrix
% R_i:      sensor noise covariance estimated in the last time instants,
%           i.e., R^i_{k/t-2}, k = t-N,...,t-1
%           USING APPROACH B: input is a single matrix
function pi_t_N = propagate_pi_i(pi_t_N_1,A_i,C_i,Q_i,R_i)
    % Retrieve dimensions
    n_i = size(A_i,1);
    o_i = size(C_i,1);
    N = size(R_i,1);
    % Compute extended observability matrix
    O_i = zeros(o_i*N,n_i);
    O_i(1:o_i,:) = C_i; 
    for k = 1:N-1
        O_i(k*o_i+1:(k+1)*o_i,:) = O_i((k-1)*o_i+1:k*o_i,:)*A_i;
    end
    % Compute blk diag R and Q
    bold_R = zeros(o_i*N);
    for k = 1:N
        %bold_R((k-1)*o_i+1:k*o_i,(k-1)*o_i+1:k*o_i) = R_i{k,1};
        bold_R((k-1)*o_i+1:k*o_i,(k-1)*o_i+1:k*o_i) = R_i;
    end
    bold_Q = zeros(n_i*(N-1));
    for k = 1:N-1
        %bold_Q((k-1)*n_i+1:k*n_i,(k-1)*n_i+1:k*n_i) = Q_i{k+1,1};
        bold_Q((k-1)*n_i+1:k*n_i,(k-1)*n_i+1:k*n_i) = Q_i;
    end
    % Compute extended controllability matrix
    C_w_N = zeros(N*o_i,(N-1)*n_i);
    for k = 2:N
        for j = 1:k-1
            C_w_N((k-1)*o_i+1:k*o_i,(j-1)*n_i+1:j*n_i) = C_i*(A_i)^(k-j-1);
        end
    end
    % Compute tilde{R}_{N/t-2}
    tilde_R = bold_R + C_w_N*bold_Q*C_w_N';    
    % Compute pi_t_N
%     pi_t_N = A_i*pi_t_N_1*A_i' + Q_i{1,1} - ...
%              A_i*pi_t_N_1*O_i'/(O_i*pi_t_N_1*O_i'+tilde_R)*O_i*pi_t_N_1*A_i';
    pi_t_N = A_i*pi_t_N_1*A_i' + Q_i - ...
             A_i*pi_t_N_1*O_i'/(O_i*pi_t_N_1*O_i'+tilde_R)*O_i*pi_t_N_1*A_i';
end

% Farina et al. PMHE1 problem 
% Inputs:
% y_i:      outputs throughout the window, i.e., y^i_k, k = t-N+1,...,t
% u^{i,y}_k:coupled output contributions from other systems, i.e.,
%           tilde{C}^[i]*hat{x}_{k/t-1}, k = t-N+1,...,t
% u^{i,x}_k:coupled dynamics contributions from other systems, i.e.,
%           tilde{A}^[i]*hat{x}_{k/t-1}, k = t-N,...,t-1
% A_i:      matrix A^[i]
% C_i:      matrix C^[i]
% Q_i:      process noise covariance estimated in the last time instants,
%           i.e., Q^i_{k/t-1}, k = t-N,...,t-1 
%           USING APPROACH B: input is a single matrix
% R_i:      sensor noise covariance estimated in the last time instants,
%           i.e., R^i_{k/t-1}, k = t-N+1,...,t
%           USING APPROACH B: input is a single matrix
% x_hat_i_t_N_t_1: hat{x}^i_{t-N/t-1}
% pi_i_t_N_t_1: pi^i_{t-N/t-1}
% Outputs:
% x:        hat{x}^i_{k/t}, k = t-N,...,t
function x = MHE_i_problem(y_i,u_i_y,u_i_x,A_i,C_i,Q_i,R_i,x_hat_i_t_N_t_1,pi_i_t_N_t_1)
    %%%% Retrieve dimensions
    N = size(y_i,2);
    n_i = size(A_i,1);
    o_i = size(R_i,1);
    %%%% Cvx for convex optimization
    cvx_begin quiet %SDPT3
    %%%% Optimization variables 
    % Initial window estimate, hat{x}^i_{t-N/t}
    variable x_hat_i_t_N(n_i)
    % Process noise throughout the window, hat{w}^i_k, k = t-N,...,t-1
    variable w_hat_i(n_i,N)
    %%%% Expressions
    % Define state estimates throught the window hat{x}^i_{k/t}, k = t-N+1,t
    expression x_hat_i_k(n_i,N)
    % First new estimate
    x_hat_i_k(:,1) = A_i*x_hat_i_t_N + u_i_x(:,1) + w_hat_i(:,1);
    % Remainder of the window
    for k = 2:N
        x_hat_i_k(:,k) = A_i*x_hat_i_k(:,k-1) + u_i_x(:,k) + w_hat_i(:,k);
    end
    % Define sensor noise throughout the window, hat{v}^i_k, k = t-N+1,...,t
    expression v_hat_i(o_i,N)
    for k = 1:N
        v_hat_i(:,k) = y_i(:,k) - C_i*x_hat_i_k(:,k) - u_i_y(:,k);
    end
    % Cost function
    expression J(1)
    % Add initial penalty
    J(1) = (1/2)*quad_form(x_hat_i_t_N-x_hat_i_t_N_t_1,inv(pi_i_t_N_t_1));
    % Add stage costs
    for k = 1:N
%         J(1) = J(1) +...
%                (1/2)*w_hat_i(:,k)'*inv(Q_i{k,1})*w_hat_i(:,k) + ... % k = t-N,...,t-1
%                (1/2)*v_hat_i(:,k)'*inv(R_i{k,1})*v_hat_i(:,k); % k = t-N+1,...,t
        J(1) = J(1) +...
               (1/2)*quad_form(w_hat_i(:,k),inv(Q_i)) + ... % k = t-N,...,t-1
               (1/2)*quad_form(v_hat_i(:,k),inv(R_i)); % k = t-N,...,t-1
    end
    %%%% Solve problem 
    minimize(J(1))
    cvx_end
    %%%% Output results
    x = [x_hat_i_t_N x_hat_i_k];
end



