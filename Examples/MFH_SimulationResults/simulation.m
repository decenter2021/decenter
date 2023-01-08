function simulation(model,methods,parameters)

    % Parameters

    % MFH_w_ss
    % PMHE1_N
    % sim_T

    %% Setup
    % For consistency
    rng(0);

    %%%% Load model 
    load(sprintf("./models/model%d/model.mat",model),'A','C','Q','R','B','N')
    n_g = size(A,1);
    o_g = size(C,1);
    m_g = size(B,2);

    %%%% Load filter synthesis
    % Store gains nicely
    K = cell(4,1);
    % Centralized
    if methods(1)
        load(sprintf("./models/model%d/synth_C.mat",model),'KC');
        K{1,1} = KC;
    end
    % OS
    if methods(2)
        load(sprintf("./models/model%d/synth_OS.mat",model),'KOS');
        K{2,1} = KOS;
    end
    % FH
    if methods(3)
        load(sprintf("./models/model%d/synth_FH.mat",model),'KFH');
        K{3,1} = KFH;
    end
    % MFH
    if methods(4)
        w_ss = parameters.MFH_w_ss;
        load(sprintf("./models/model%d/synth_MFH_%d.mat",model,w_ss),'KMFH');
        K{4,1} = KMFH;
    end

    %%%% Load controller (it is introduced just to avoid numerical problems)
    load(sprintf("./models/model%d/synth_ctrl.mat",model),'KCtrlOS');

    % Simulation parameters
    SimIt = parameters.sim_T;
    
    %% Ground truth simulation
    % I could simulate just the error dynamics for C, OS, FH, and MFH
    % But I have to simulate the state evolution for PMHE1, so I do it like
    % that to make it easier to have equal process and sensor noise
    %%%% Init simulation variables 
    x = zeros(n_g,SimIt+1); % x_t, t = 0,...,SimIt
    u = zeros(m_g,SimIt); % u_t, t = 0,...,SimIt-1
    y = zeros(o_g,SimIt); % y_t, t = 1,...,SimIt
    %%%% Initial simulation parameters
    Px0 = 10*eye(n_g);
    x(:,1) = mvnrnd(zeros(n_g,1),Px0)'; 
    %%%% Ground-truth simulation
    for k = 1:SimIt
        u(:,k) = -KCtrlOS*x(:,k);
        x(:,k+1) = A*x(:,k) + B*u(:,k) + mvnrnd(zeros(n_g,1),Q)';
        % Warning: the indices of x and y start at different time-instants
        y(:,k) = C*x(:,k+1) + mvnrnd(zeros(o_g,1),R)'; 
    end

    %% Filter simulation
    %%%% Init filter simulation variables
    x_hat = cell(5,1);

    %%%% Filter simulations initialization
    P0 = 2*eye(size(A)); %zeros(n_g,n_g);
    x_hat_0 = x(:,1) + mvnrnd(zeros(n_g,1),P0)';
   
    %%%% Luenberger filters: C, OS, FH
    for m = 1:3
        if ~methods(m), break; end
        % Init filter simulation matrices
        x_hat{m,1} = zeros(n_g,SimIt+1);
        % Initial estimate 
        x_hat{m,1}(:,1) = x_hat_0;
        % Simulate filter dynamics
        for k = 1:SimIt
            % Predict: Compute x_{t/t-1}
            x_aux = A*x_hat{m,1}(:,k) + B*u(:,k);
            % Filter with y_t: Compute x_{t/t}
            x_hat{m,1}(:,k+1) = x_aux + K{m,1}*(y(:,k)-C*x_aux);
        end
    end
    fprintf("Luenberger filters simulation completed.\n");

    %%%% MFH 
    m = 4; % MFH
    if methods(m)
        KMFH = K{m,1};
        % Init filter simulation matrices
        x_hat{m,1} = zeros(n_g,SimIt+1);
        % Initial estimate 
        x_hat{m,1}(:,1) = x_hat_0;
        % Simulate filter dynamics
        for k = 1:SimIt
            % If it is not possible to use full steady-state window, use OS
            if k < w_ss
                % Predict: Compute x_{t/t-1}
                x_aux = A*x_hat{m,1}(:,k) + B*u(:,k);
                % Filter with y_t: Compute x_{t/t}
                x_hat{m,1}(:,k+1) = x_aux + K{2,1}*(y(:,k)-C*x_aux);
            else
                x_aux = x_hat{m,1}(:,k-w_ss+1);
                for l = k-w_ss+1:k
                    x_aux = A*x_aux + B*u(:,l);
                    x_aux = x_aux + KMFH{l-(k-w_ss),1}*(y(:,l)-C*x_aux);
                end
                x_hat{m,1}(:,k+1) = x_aux;
            end
        end
        fprintf("MFH filter simulation completed.\n");
    end


    %%%% PMHE1 - Farina et al. 2010
    m = 5;
    if methods(m)
        % PMHE1 parameters 
        PMHE1_N = parameters.PMHE1_N;
        % Retrieve dimensions
        n_i = size(A,1)/N;
        % Init filter simulation matrices
        x_hat{m,1} = zeros(n_g,SimIt+1);
        % Initial estimate 
        x_hat{m,1}(:,1) = x_hat_0;
        % Initial pi_t
        pi_i_t_N_1_t_2 = cell(N,1);
        for i = 1:N
            pi_i_t_N_1_t_2{i,1} = P0((i-1)*n_i+1:i*n_i,(i-1)*n_i+1:i*n_i);
        end
        % Simulation 
        for k = 1:SimIt
            fprintf("PMHE1 filter iteration: %d\n",k);
            % x_hat_t_1: hat{x}_{k/t-1}, k = t-N,...,t-1
            % y_PMHE1: y^i_k, k = t-N+1,...,t
            % u_PMHE1: u_k, k = t-N,...,t-1  
            if k == 1
                % Window of length 1
                x_hat_t_1 = x_hat_0;
                y_PMHE1 = y(:,1);
                u_PMHE1 = u(:,1);
            elseif k <= PMHE1_N
                % Window of length k
                x_hat_t_1 = x_hat_i_t;
                y_PMHE1 = y(:,1:k);
                u_PMHE1 = u(:,1:k);
            else
                x_hat_t_1 = x_hat_i_t(:,2:end);
                y_PMHE1 = y(:,k-PMHE1_N+1:k);
                u_PMHE1 = u(:,k-PMHE1_N+1:k);
            end
            [x_hat_i_t,pi_i_t_N_1_t_2] = farina_et_al_2010(A,C,Q,R,B,x_hat_t_1,y_PMHE1,u_PMHE1,pi_i_t_N_1_t_2);
            % x: hat{x}^i_{k/t}, k = t-N,...,t
            x_hat{m,1}(:,k+1) = x_hat_i_t(:,end);
        end
        fprintf("PMHE1 filter simulation completed.\n");
    end

    %% Post processing
    error = cell(5,1);
    error_norm = zeros(5,SimIt);
    for m = 1:5
        if ~methods(m), continue; end
        error{m,1} = x_hat{m,1}-x;
        for k = 1:SimIt+1
            error_norm(m,k) = norm(error{m,1}(:,k));
        end
    end
    save(sprintf("./models/model%d/simulation.mat",model),'error','x','x_hat','methods','parameters','error_norm');
    
end