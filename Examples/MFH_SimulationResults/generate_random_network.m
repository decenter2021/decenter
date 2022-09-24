function generate_random_network(seed)

    %% Model parameters
    % Graph generation parameters 
    % ---------------------------
    rng(seed);
    N = 20;
    max_neighborhood = 5;
    proximity_neighborhood = 7;
    q_rel = 0.5;
    r_rel = 0.5;
    a_rel_out_diag = 0.5;
    a_rel_coupling = 0.5;
    q_rel_coupling = 0.5;
    
    % ---------------------------
    % Model 1:
    % rng(seed);
    % N = 20;
    % max_neighborhood = 5;
    % proximity_neighborhood = 7;
    % q_rel = 0.5;
    % r_rel = 0.5;
    % a_rel_out_diag = 0.5;
    % a_rel_coupling = 0.5;
    % q_rel_coupling = 0.5;

    % Model 40:
    % rng(seed);
    % N = 500;
    % max_neighborhood = 5;
    % proximity_neighborhood = 7;
    % q_rel = 0.5;
    % r_rel = 0.5;
    % a_rel_out_diag = 0.3;
    % a_rel_coupling = 0.4;
    % q_rel_coupling = 0.3;

    % Model 41:
    % rng(40);
    % N = 500;
    % max_neighborhood = 5;
    % proximity_neighborhood = 7;
    % q_rel = 0.5;
    % r_rel = 0.5;
    % a_rel_out_diag = 0.3;
    % a_rel_coupling = 0.4;
    % q_rel_coupling = 0;

    % Model 42:
    % rng(seed);
    % N = 1000;
    % max_neighborhood = 5;
    % proximity_neighborhood = 7;
    % q_rel = 0.5;
    % r_rel = 0.5;
    % a_rel_out_diag = 0.4;
    % a_rel_coupling = 0.05;
    % q_rel_coupling = 0;

    % Model 47
    % N = 500;
    % max_neighborhood = 5;
    % proximity_neighborhood = 7;
    % q_rel = 0.5;
    % r_rel = 0.5;
    % a_rel_out_diag = 0.4;
    % a_rel_coupling = 0.4;
    % q_rel_coupling = 0.4;

    %% Generate model
    edges = [];    
    % Global dynamics
    n = 2;
    o = 1;
    m = 1;
    A = zeros(N*n);
    C = zeros(o*N,n*N);
    Q = zeros(N*n);
    R = zeros(o*N,o*N);
    B = zeros(n*N,m*N);
    itGeneration = 0;
    
    while true
        itGeneration = itGeneration +1;
        for i = 1:N
            % Principal diagonal A,Q
            A((i-1)*n+1:(i-1)*n+n,(i-1)*n+1:(i-1)*n+n) = 0.8*eye(n)+a_rel_out_diag*(rand(n,n)-0.5);
            Q((i-1)*n+1:(i-1)*n+n,(i-1)*n+1:(i-1)*n+n) = ...
                q_rel*(A((i-1)*n+1:(i-1)*n+n,(i-1)*n+1:(i-1)*n+n)*...
                A((i-1)*n+1:(i-1)*n+n,(i-1)*n+1:(i-1)*n+n)');
            % Diagonal 1 A,Q
            lcur = [];
            lmax =  round(exp(rand()*(log(max_neighborhood))));
            l = 0;
            while true
                if l == lmax
                    break;
                end
                l = l+1;
                j = round(rand()*(N-1))+1;
                if j ~= i 
                    if ~sum(lcur == j) && abs(j-i) <= proximity_neighborhood 
                        A((i-1)*n+1:(i-1)*n+n,(j-1)*n+1:(j-1)*n+n) = a_rel_coupling*(rand(n,n)-0.5);
                        Q((i-1)*n+1:(i-1)*n+n,(j-1)*n+1:(j-1)*n+n) = ...
                            q_rel_coupling*(A((i-1)*n+1:(i-1)*n+n,(j-1)*n+1:(j-1)*n+n)*...
                            A((j-1)*n+1:(j-1)*n+n,(i-1)*n+1:(i-1)*n+n)');
                        Q((j-1)*n+1:(j-1)*n+n,(i-1)*n+1:(i-1)*n+n) = ...
                           Q((i-1)*n+1:(i-1)*n+n,(j-1)*n+1:(j-1)*n+n)';
                        edges = [edges [j;i]];
                        lcur = [lcur j];
                    else
                        l =l-1;
                        continue;
                    end
                else
                    l = l-1;
                    continue;
                end
    
            end
    
            % C
            C((i-1)*o+1:(i-1)*o+o,(i-1)*n+1:(i-1)*n+n) = 2*(rand(o,n)-0.5);
            % R 
            R((i-1)*o+1:(i-1)*o+o,(i-1)*o+1:(i-1)*o+o) = ...
                r_rel*(C((i-1)*o+1:(i-1)*o+o,(i-1)*n+1:(i-1)*n+n)*...
                C((i-1)*o+1:(i-1)*o+o,(i-1)*n+1:(i-1)*n+n)');
        end
        if min(eig(Q)) > 0
            break;
        elseif itGeneration == 100
            fprintf("Exceed the maximum number of iterations\n"); 
            exit();
        end
    end
    
    % Decoupled matrix B
    for i = 1:N
        B((i-1)*n+1:i*n,(i-1)*m+1:i*m) = rand(n,m);
    end

    % Sparsity pattern 
    E = zeros(n*N,o*N);
    for i = 1:N
        % C
        E((i-1)*n+1:(i-1)*n+n,(i-1)*o+1:(i-1)*o+o) = ones(n,o);
    end
    
    %% Plot diagraph
    G = digraph();
    for e = 1:size(edges,2)
        G = addedge(G,edges(1,e),edges(2,e));
    end
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    set(gca,'FontSize',35);
    ax = gca;
    axis off;
    p = plot(G,'Layout','force','UseGravity',true,'Iterations',1e3);
    p.Marker = 's';
    p.NodeColor = 'r';
    p.MarkerSize = 20; %10 %20;
    p.LineWidth = 2;
    p.ArrowSize = 20; %10
    p.NodeLabel = {};
    hold off;
    mkdir(sprintf("./models/model%d",seed));
    diagraph_path = sprintf("./models/model%d/diagraph.fig",seed);
    savefig(diagraph_path)
    close all;

    %% Save model  
    model_path = sprintf("./models/model%d/model.mat",seed);
    save(model_path);

    %% Model properties
    fprintf("Global state size: %d\n", n*N)
    fprintf("Maximum absolute eigenvalue: %g\n", max(abs(eig(A))));

end
