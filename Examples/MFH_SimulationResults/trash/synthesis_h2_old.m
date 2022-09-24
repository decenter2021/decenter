% Makes use of:
% YALMIP [1] (https://yalmip.github.io) 
% SEDUMI [2] (https://github.com/sqlp/sedumi)

function [K,P,mu] = synthesis_h2(seed)
    % Load model 
    model_path = sprintf("./models/model%d/model.mat",seed);
    load(model_path,'A','C','Q','R','E')
    % Load centralized performance
    load(sprintf("./models/model%d/synth_C.mat",seed),'PC');
    % Precompute sqr of process and sensor noise
    Q1_2 = sqrtm(Q);
    R1_2 = sqrtm(R);
    % Load OS initialization
    load(sprintf("./models/model%d/synth_OS.mat",seed),'KOS');
    % Init P-K iterations
    muPrev = inf;
    K = (rand(size(E))-0.5).*E;
    fprintf("Initialization: OS method.\n")
    fprintf("Maximum absolute CL eigenvalue: %f\n",max(abs(eig((A-K*C)))));
    countMax = 100;
    epsl_inf = 1e-3;
    % P-K iterations
    count = 0;
    tic;
    while true
        count = count +1;
        fprintf("Iteration %d.\n",count);
        if count > countMax
           mu = nan;
           break;
        end
        [~,P] = optimizeP(A,C,Q1_2,R1_2,K);
        [mu,K] = optimizeK(A,C,Q1_2,R1_2,E,P);
        fprintf("Mu: %g.\n",mu);
        fprintf("Maximum absolute CL eigenvalue: %f\n",max(abs(eig((A-K*C)))));
        if abs(mu-muPrev)/muPrev < epsl_inf
            break;
        end
        muPrev = mu;
    end
    toc;
    if ~isnan(K)
        fprintf("Trace H2: %f\n",trace(P));
        fprintf("Trace H2/C: %f\n",trace(P)/trace(PC));
        fprintf("Maximum absolute CL eigenvalue: %f\n",max(abs(eig((A-K*C)))));
        save(sprintf("./models/model%d/synth_H2.mat",seed),'K','P','mu');
    end
end

function [mu,Popt] = optimizeP(A,C,Q1_2,R1_2,K)
    % Get size from variables
    n = size(A,1);
    o = size(C,1);
    % Optimization variables
    P = sdpvar(n); 
    W = sdpvar(n); 
    mu = sdpvar(1);
    % LMI matrices
    X = [P                      (A-K*C)*P       [Q1_2 -K*R1_2];...
         P*(A-K*C)'   P               zeros(n,n+o);...
         [Q1_2 -K*R1_2]'        zeros(n,n+o)'   eye(n+o)];
    W_ = [W P; P P];
    % Contraints 
    epsl = 1e-5;
    constraints = [ P >= epsl*eye(size(P,1));...
                    W >= epsl*eye(size(W,1));...
                    mu >= epsl*eye(size(mu,1));...
                    X >= epsl*eye(size(X,1));...
                    W_ >= epsl*eye(size(W_,1));...
                    trace(W) <= mu;];
    % Solve
    ops = sdpsettings('verbose',0,'warning',1,'solver','sedumi');
    optimize(constraints,mu,ops);
    Popt = value(P);
    mu = value(mu);
end

function [mu,Kopt] = optimizeK(A,C,Q1_2,R1_2,E,P)
    % Get size from variables
    n = size(A,1);
    o = size(C,1);
    % Optimization variables
    W = sdpvar(n); 
    K = sdpvar(n,o,'full');
    mu = sdpvar(1);
    % LMI matrices
    X = [P                      (A-K*C)*P       [Q1_2 -K*R1_2];...
         P*(A-K*C)'   P               zeros(n,n+o);...
         [Q1_2 -K*R1_2]'        zeros(n,n+o)'   eye(n+o)];
    W_ = [W P; P P];
    % Contraints 
    epsl = 1e-5;
           
    constraints = [ K(E==0) == 0;...
                    W >= epsl*eye(size(W,1));...
                    mu >= epsl*eye(size(mu,1));...
                    X >= epsl*eye(size(X,1));...
                    W_ >= epsl*eye(size(W_,1));...
                    trace(W) <= mu;];
    % Solve 
    ops = sdpsettings('verbose',0,'warning',1,'solver','sedumi');
    optimize(constraints,mu,ops);
    Kopt = value(K);
    mu = value(mu);
end

% References 
% [1] Lofberg, J., 2004, September. YALMIP: A toolbox for modeling and 
% optimization in MATLAB. In 2004 IEEE international conference on robotics
% and automation (IEEE Cat. No. 04CH37508) (pp. 284-289). IEEE.
% [2] Sturm, J.F., 1999. Using SeDuMi 1.02, a MATLAB toolbox for 
% optimization over symmetric cones. Optimization methods and software, 
% 11(1-4), pp.625-653.


