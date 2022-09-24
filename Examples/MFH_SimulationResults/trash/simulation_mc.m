function simulation_mc(model,methods,parameters)

    % Parameters
    % MFH_w_ss
    % sim_T
    % PMHE1_N
    % mc_it
    
    % MFH method is initialized with OS methods
    if methods(5), methods(2) = true; end

    %% Setup
    %%%% Load model 
    load(sprintf("./models/model%d/model.mat",model),'A','C','Q','R','B','N');
    n_g = size(A,1);
    o_g = size(C,1);
    m_g = size(B,2);

    %%%% Load filter synthesis
    % Store gains nicely
    K = cell(5,1);
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
    % H2 
    if methods(4)
        w_ss = parameters.MFH_w_ss;
        data = load(sprintf("./models/model%d/synth_H2.mat",model),'K');
        K{4,1} = data.K;
    end

    % MFH
    if methods(5)
        w_ss = parameters.MFH_w_ss;
        load(sprintf("./models/model%d/synth_MFH_%d.mat",model,w_ss),'KMFH');
        K{5,1} = KMFH;
    end

    %%%% Load controller (it is introduced just to avoid numerical problems)
    load(sprintf("./models/model%d/synth_ctrl.mat",model),'KCtrlOS');

    % Simulation parameters
    SimIt = parameters.sim_T;

    %% Compute theoretical performance throughout simulation
    
    % Init evolution of P
    P_ev = cell(5,SimIt);
    P_ss = cell(5,1);
    P_tr_ss_epsl = zeros(4,1);
    P0 = zeros(size(A));

    % Cent, OS, FH, H2
    for m = 1:4
        if ~methods(m), continue; end
        for i = 1:SimIt
            if i > 1 
                P_ = A*P_ev{m,i-1}*transpose(A)+Q;
            else
                P_ = A*P0*transpose(A)+Q;
            end
            P_ev{m,i} = K{m,1}*R*K{m,1}'+(eye(n_g)-K{m,1}*C)*P_*(eye(n_g)-K{m,1}*C)';
        end
        P_ss{m,1} = P_ev{m,end};
        P_tr_ss_epsl(m) = (trace(P_ev{m,end})-trace(trace(P_ev{m,end-1})))/trace(P_ev{m,end});
    end
    m = 5; % MFH
    if methods(m)
        MFH_Wss = parameters.MFH_w_ss;
        KMFH = K{m,1};
        for i = 1:SimIt
            if i < MFH_Wss
                P_ev{m,i} = P_ev{2,i}; % MFH initialized with one-step
            else
                if i ~= MFH_Wss
                    P_ = A*P_ev{m,i-MFH_Wss}*transpose(A)+Q;
                else
                    P_ = A*P0*transpose(A)+Q;
                end
                for j = i-MFH_Wss+1:i
                    Paux = KMFH{j-i+MFH_Wss,1}*R*transpose(KMFH{j-i+MFH_Wss,1})+...
                        (eye(n_g)-KMFH{j-i+MFH_Wss,1}*C)*P_*...
                        transpose(eye(n_g)-KMFH{j-i+MFH_Wss,1}*C);
                    P_ = A*Paux*transpose(A)+Q;
                end
                P_ev{m,i} = Paux;
            end
        end
        P_ss{m,1} = P_ev{m,end};
        P_tr_ss_epsl(m) = (trace(P_ev{m,end})-trace(trace(P_ev{m,end-1})))/trace(P_ev{m,end});
    end
    
    if methods(1)
        fprintf("Trace C: \t%f \t(ri: %g)\n",trace(P_ss{1,1}),P_tr_ss_epsl(1));
        if methods(2)
            fprintf("Trace OS/C: \t%f \t(ri: %g)\n",trace(P_ss{2,1})/trace(P_ss{1,1}),P_tr_ss_epsl(2));
        end
        if methods(3)
            fprintf("Trace FH/C: \t%f \t(ri: %g)\n",trace(P_ss{3,1})/trace(P_ss{1,1}),P_tr_ss_epsl(3));
        end
        if methods(4)
            fprintf("Trace H2/C: \t%f\t(ri: %g)\n",trace(P_ss{4,1})/trace(P_ss{1,1}),P_tr_ss_epsl(4));
        end
        if methods(5)
            fprintf("Trace MFH/C: \t%f\t(ri: %g)\n",trace(P_ss{5,1})/trace(P_ss{1,1}),P_tr_ss_epsl(5));
        end
    else      
        if methods(2)
            fprintf("Trace OS: \t%f \t(ri: %g)\n",trace(P_ss{2,1}),P_tr_ss_epsl(2));
        end
        if methods(3)
            fprintf("Trace FH: \t%f \t(ri: %g)\n",trace(P_ss{3,1}),P_tr_ss_epsl(3));
        end
        if methods(4)
            fprintf("Trace H2: \t%f \t(ri: %g)\n",trace(P_ss{4,1}),P_tr_ss_epsl(4));
        end
        if methods(5)
            fprintf("Trace MFH: \t%f \t(ri: %g)\n",trace(P_ss{5,1}),P_tr_ss_epsl(5));
        end
    end

    % Save theoretical data
    save(sprintf("./models/model%d/projected_ss.mat",model),'P_ss','methods','parameters');

    %% MC simulations

    error = cell(6,1);
    for m = 1:6
        if ~methods(m), continue; end
        % Terrible code, I know, but it runs much faster :)
        if m == 1
            errorC = zeros(n_g,SimIt+1,parameters.mc_it);
        elseif m == 2
            errorOS = zeros(n_g,SimIt+1,parameters.mc_it);
        elseif m == 3
            errorFH = zeros(n_g,SimIt+1,parameters.mc_it);
        elseif m == 4
            errorH2 = zeros(n_g,SimIt+1,parameters.mc_it);
        elseif m == 5
            errorMFH = zeros(n_g,SimIt+1,parameters.mc_it);
        else
            errorPMHE1 = zeros(n_g,SimIt+1,parameters.mc_it);
        end
    end

    for i = 1:parameters.mc_it
    %parfor i = 1:parameters.mc_it
        fprintf("Monte Carlo simulation %d/%d.\n",i,parameters.mc_it);
        error_it = simulation_mc_it(model,methods,parameters,A,C,Q,R,B,N,K,KCtrlOS);
        for m = 1:6
            if ~methods(m), continue; end
            % Terrible code, I know, but it runs much faster :)
            if m == 1
                errorC(:,:,i) = error_it{m,1};
            elseif m == 2
                errorOS(:,:,i) = error_it{m,1};
            elseif m == 3
                errorFH(:,:,i) = error_it{m,1};
            elseif m == 4   
                errorH2(:,:,i) = error_it{m,1};
            elseif m == 5
                errorMFH(:,:,i) = error_it{m,1};
            else
                errorPMHE1(:,:,i) = error_it{m,1};
            end
        end
    end

    fprintf("Finished Monte Carlo simulations.\n");

    for m = 1:6
        if ~methods(m), continue; end
        % Terrible code, I know, but it runs much faster :)
        if m == 1
            error{m,1} = errorC;
            clear errorC;
        elseif m == 2
            error{m,1} = errorOS;
            clear errorOS;
        elseif m == 3
            error{m,1} = errorFH;
            clear errorFH;
        elseif m == 4
            error{m,1} = errorH2;            
            clear errorMFH;
        elseif m == 5
            error{m,1} = errorMFH;
            clear errorMFH;
        else
            error{m,1} = errorPMHE1;
            clear errorPMHE1;
        end
    end

    save(sprintf("./models/model%d/mc_raw_error.mat",model),'error','methods','parameters');

    fprintf("Processing simulations.\n");
    P_mc = cell(6,SimIt+1);
    for m = 1:6       
        if ~methods(m), continue; end
        fprintf("Processing method %d/%d.\n",m,6);
        %parfor t = 1:SimIt+1
        for t = 1:SimIt+1          
            P_mc{m,t} = cov(permute(error{m,1}(:,t,:),[1 3 2])');
        end
    end

    clear error;

    tr_P_mc = zeros(6,SimIt+1);
    for m = 1:6
        if ~methods(m), continue; end
        for t = 1:SimIt+1
            tr_P_mc(m,t) = trace(P_mc{m,t});
        end
    end

    save(sprintf("./models/model%d/mc.mat",model),'P_mc','tr_P_mc','methods','parameters');
    

end


function error = simulation_mc_it(model,methods,parameters,A,C,Q,R,B,N,K,KCtrlOS)

    % Parameters
    % MFH_w_ss
    % PMHE1_N
    % sim_T

    n_g = size(A,1);
    o_g = size(C,1);
    m_g = size(B,2);

    % Simulation parameters
    SimIt = parameters.sim_T;
    w_ss = parameters.MFH_w_ss;
    
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
    x_hat = cell(6,1);

    %%%% Filter simulations initialization
    P0 = 0*eye(size(A)); %zeros(n_g,n_g);
    x_hat_0 = x(:,1) + mvnrnd(zeros(n_g,1),P0)';
   
    %%%% Luenberger filters: C, OS, FH
    for m = 1:4
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

    %%%% MFH 
    m = 5; % MFH
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
    end


    %%%% PMHE1 - Farina et al. 2010
    m = 6;
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
    end

    %% Post processing
    error = cell(6,1);
    for m = 1:6
        if ~methods(m), continue; end
        error{m,1} = x_hat{m,1}-x;
    end
end




