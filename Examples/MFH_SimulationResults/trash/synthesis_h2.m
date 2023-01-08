% Makes use of:
% YALMIP [1] (https://yalmip.github.io) 
% SEDUMI [2] (https://github.com/sqlp/sedumi)

function [K,P,mu] = synthesis_h2(seed)
    % Load model 
    model_path = sprintf("./models/model%d/model.mat",seed);
    load(model_path,'A','C','Q','R','E')
    % Load centralized performance
    load(sprintf("./models/model%d/synth_C.mat",seed),'PC','KC');
    % Precompute sqr of process and sensor noise
    Q1_2 = sqrtm(Q);
    R1_2 = sqrtm(R);
%     % Load OS initialization
%     load(sprintf("./models/model%d/synth_OS.mat",seed),'KOS');
%     % Init P-K iterations
%     muPrev = inf;
%     K = KC.*E;
%     fprintf("Initialization: random gain.\n")
%     fprintf("Maximum absolute CL eigenvalue: %f\n",max(abs(eig((A-K*C)))));
%     countMax = 100;
%     epsl_inf = 1e-3;
%     % P-K iterations
%     count = 0;
    tic;
%     while true
%         count = count +1;
%         fprintf("Iteration %d.\n",count);
%         if count > countMax
%            mu = nan;
%            break;
%         end
%         P = optimizeP(A,C,Q1_2,R1_2,K);
%         [mu,K] = optimizeK(A,C,Q1_2,R1_2,E,P,K);
%         fprintf("Mu: %g.\n",mu);
%         fprintf("Maximum absolute CL eigenvalue: %f\n",max(abs(eig((A-K*C)))));
%         if abs(mu-muPrev)/muPrev < epsl_inf
%             break;
%         end
%         muPrev = mu;
%     end   
    [mu,K,P] = optimizeAll(A,C,Q1_2,R1_2,E);
    toc;
    if ~isnan(K)
        fprintf("Trace H2: %f\n",trace(P));
        fprintf("Trace H2/C: %f\n",trace(P)/trace(PC));
        fprintf("Maximum absolute CL eigenvalue: %f\n",max(abs(eig((A-K*C)))));
        save(sprintf("./models/model%d/synth_H2.mat",seed),'K','P','mu');
    end
end

function Popt = optimizeP(A,C,Q1_2,R1_2,K)
    % Get size from variables
    n = size(A,1);
    o = size(C,1);
    % Optimization variables
    P = sdpvar(n); 
    W = sdpvar(n); 
    mu = sdpvar(1);
    % LMI matrices
    X = [P                              (eye(n)-K*C)*A*P    [(eye(n)-K*C)*Q1_2 K*R1_2];...
         ((eye(n)-K*C)*A*P)'            P                   zeros(n,n+o);...
         [(eye(n)-K*C)*Q1_2 K*R1_2]'    zeros(n,n+o)'       eye(n+o)];
    W_ = [W P; P' P];
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
end

function [mu,Kopt] = optimizeK(A,C,Q1_2,R1_2,E,P,K1)
    % Get size from variables
    n = size(A,1);
    o = size(C,1);
    % Optimization variables
    W = sdpvar(n); 
    %[i,j,~] = find(E);
    %K = sparse(i,j,sdpvar(length(i),1));
    K = sdpvar(n,o,'full');
    mu = sdpvar(1);
    % LMI matrices
    X = [P                              (eye(n)-K*C)*A*P    [(eye(n)-K1*C)*Q1_2 K1*R1_2];...
         ((eye(n)-K*C)*A*P)'            P                   zeros(n,n+o);...
         [(eye(n)-K1*C)*Q1_2 K1*R1_2]'    zeros(n,n+o)'       eye(n+o)];
    W_ = [W P; P' P];
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

function [mu,Kopt,Popt] = optimizeAll(A,C,Q1_2,R1_2,E)
    % Get size from variables
    n = size(A,1);
    o = size(C,1);

    %define decision variables
    P = sdpvar(n);
    Z = sdpvar(n);
    Gd = sdpvar(n,o,'full');
    v = sdpvar(1);

    Ad = A;
    Cd1 = eye(n);
    Cd2 = C;
    Bd1 = [eye(n) zeros(n,o)];
    Dd21 = [zeros(o,n) eye(o)];

    %Build LMI
    Mat1 = [P (P*Ad - Gd*Cd2) (P*Bd1 - Gd*Dd21);
        (P*Ad - Gd*Cd2)' P zeros(size((P*Bd1 - Gd*Dd21)));
        (P*Bd1 - Gd*Dd21)' zeros(size((P*Bd1 - Gd*Dd21)))' eye(n+o)];
    Mat2 = [Z P*Cd1;
        (P*Cd1)' P];
    K = inv(P)*Gd;

%define constraints
F = [];
F = [F; K(E==0) == 0];
F = [F; Mat1>=smallnum*eye(size(Mat1))];
F = [F; Mat2>=smallnum*eye(size(Mat2))];
F = [F; trace(Z)<=v];
F = [F; P >=smallnum*eye(size(P))];
F = [F; Z >=smallnum*eye(size(Z))];

%optimize
solution = optimize(F,v);

%print optimal H2 observer
Kopt = inv(value(P))*value(Gd);
Popt = value(P);
%print out H2 norm
mu = sqrt(value(v));

%     % Optimization variables
%     W = sdpvar(n); 
%     P = sdpvar(n); 
%     [i,j,~] = find(E);
%     K = sparse(i,j,sdpvar(length(i),1));
%     %K = sdpvar(n,o,'full');
%     mu = sdpvar(1);
%     % LMI matrices
%     X = [P                              (eye(n)-K*C)*A*P    [(eye(n)-K*C)*Q1_2 K*R1_2];...
%          ((eye(n)-K*C)*A*P)'            P                   zeros(n,n+o);...
%          [(eye(n)-K*C)*Q1_2 K*R1_2]'    zeros(n,n+o)'       eye(n+o)];
%     W_ = [W P; P' P];
%     % Contraints 
%     epsl = 1e-5;
%     constraints = [ P >= epsl*eye(size(P,1));...
%                     W >= epsl*eye(size(W,1));...
%                     mu >= epsl*eye(size(mu,1));...
%                     X >= epsl*eye(size(X,1));...
%                     W_ >= epsl*eye(size(W_,1));...
%                     trace(W) <= mu;];
%     % Solve 
%     %ops = sdpsettings('verbose',1,'warning',1,'solver','sedumi');
%     optimize(constraints,mu);
%     Kopt = value(K);
%     mu = value(mu);
%     Popt = value(P);
end

% References 
% [1] Lofberg, J., 2004, September. YALMIP: A toolbox for modeling and 
% optimization in MATLAB. In 2004 IEEE international conference on robotics
% and automation (IEEE Cat. No. 04CH37508) (pp. 284-289). IEEE.
% [2] Sturm, J.F., 1999. Using SeDuMi 1.02, a MATLAB toolbox for 
% optimization over symmetric cones. Optimization methods and software, 
% 11(1-4), pp.625-653.


