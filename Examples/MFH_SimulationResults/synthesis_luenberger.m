function synthesis_luenberger(seed)
    % Load model 
    model_path = sprintf("./models/model%d/model.mat",seed);
    load(model_path,'A','C','Q','R','E','B','N')
    
    % Centralized
    optsC = struct('verbose',true,'epsl',1e-5,'maxIt',1e3);
    tic;
    [KC,PC] = kalmanCentralizedLTI(A,C,Q,R,optsC);
    toc;
    fprintf("Trace centralized: %f\n",trace(PC));
    save(sprintf("./models/model%d/synth_C.mat",seed),'KC','PC','optsC')
    
    % One-step
    optsOS = struct('verbose',true,'epsl',1e-4,'maxIt',1e4);
    tic;
    [KOS,POS] = kalmanOneStepLTI(A,C,Q,R,E,optsOS);
    toc;
    fprintf("Trace OS: %f\n",trace(POS));
    fprintf("Trace OS/C: %f\n",trace(POS)/trace(PC));
    save(sprintf("./models/model%d/synth_OS.mat",seed),'KOS','POS','optsOS')

    % Sub-optimal decentralized controller just to prevent numerical problems with
    % unstable systems
    optsCtrlOS = struct('verbose',true,'epsl',1e-4,'maxIt',1e3);
    tic;
    [KCtrlOS,PCtrlOS] = LQROneStepLTI(A,B,eye(size(A,1)),100*eye(size(B,2)),kron(eye(N),ones(size(B,2)/N,size(A,1)/N)),optsCtrlOS);
    toc;
    fprintf("Maximum absolute control CL eigenvalue: %f\n",max(abs(eig(A-B*KCtrlOS))));
    save(sprintf("./models/model%d/synth_ctrl.mat",seed),'KCtrlOS','PCtrlOS','optsCtrlOS')


end