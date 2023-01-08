function flag = farina_et_al_2010_convergence(seed,N)

    %%%% Load model 
    load(sprintf("./models/model%d/model.mat",seed),'A','C','Q','R','B','n','o')
    n_g = size(A,1);
    o_g = size(C,1);
    m_g = size(B,2);
    n_i = n;
    o_i = o;
    M = size(A,1)/n_i;

    %%%% Compute A^{star}, tilde{A}, C^{star}, tilde{C}
    %%%% And check Assumption 1
    A_star = zeros(n_g,n_g);
    C_star = zeros(o_g,n_g);
    for i = 1:M
        A_i = A((i-1)*n_i+1:i*n_i,(i-1)*n_i+1:i*n_i);
        C_i = C((i-1)*o_i+1:i*o_i,(i-1)*n_i+1:i*n_i);
        A_star((i-1)*n_i+1:i*n_i,(i-1)*n_i+1:i*n_i) = A_i;
        C_star((i-1)*o_i+1:i*o_i,(i-1)*n_i+1:i*n_i) = C_i;
        % Check Assumption 1
        if rank(obsv(A_i,C_i)) < n_i
            fprintf("Subsystem %d does not satisfy Assumption 1.\n",i);
            flag = 0;
            return;
        end
    end
    A_tilde = A-A_star;
    C_tilde = C-C_star;

    
    %%%% Compute C_{N+1}
    C_N_1 = zeros(o_g*(N+1),n_g*(N+1));
    for i = 1:N+1
        for j = 1:i
            if i == j
                C_N_1((i-1)*o_g+1:i*o_g,(j-1)*n_g+1:j*n_g) = C_tilde;
            else
                C_N_1((i-1)*o_g+1:i*o_g,(j-1)*n_g+1:j*n_g) = C_star*A_star^(i-j-1)*A_tilde;
            end
            
        end
    end

    %%%% Compute O^{star}_{N+1}
    O_star_N_1 = zeros(M*o_i*(N+1),M*n_i);
    for i = 1:M
        O_i = zeros(o_i*(N+1),n_i);
        A_i = A((i-1)*n_i+1:i*n_i,(i-1)*n_i+1:i*n_i);
        C_i = C((i-1)*o_i+1:i*o_i,(i-1)*n_i+1:i*n_i);
        for k = 1:N+1
            O_i((k-1)*o_i+1:o_i*k,:) = C_i*A_i^(k-1);
        end
        O_star_N_1((i-1)*o_i*(N+1)+1:i*o_i*(N+1),(i-1)*n_i+1:i*n_i) = O_i;
    end
    

    %%%% Compute M_1, M_2
    M1 = zeros(n_g*(N+1),n_g);
    for k = 1:N+1
        M1((k-1)*n_g+1:n_g*k,:) = A_star^k;
    end
    M2 = zeros(n_g*(N+1),n_g*(N+1));
    for i = 1:N+1
        for j = 1:i
            M2((i-1)*n_g+1:n_g*i,(j-1)*n_g+1:n_g*j) = A_star^(i-j)*A_tilde;
        end
    end

    %%%% Compute Phi_1 and check whether it is Shur
    Phi_1 = M2-M1/(O_star_N_1'*O_star_N_1)*O_star_N_1'*C_N_1;
    if max(abs(eig(Phi_1))) >= 1
        flag = 0;
        fprintf("PMHE1 is not guaranteed to converge.\n");
    else
        flag = 1;
        fprintf("PMHE1 is guaranteed to converge.\n")
    end

end
