function synthesis_fh(seed,epsl,W,maxOLIt)
    % Load model 
    model_path = sprintf("./models/model%d/model.mat",seed);
    load(model_path,'A','C','Q','R','E')

    % Load centralized performance
    load(sprintf("./models/model%d/synth_C.mat",seed),'PC');
    % Load LTI OS gain
    load(sprintf("./models/model%d/synth_OS.mat",seed),'POS');

    % Finite-horizon
    optsFH = struct('verbose',true,'epsl',epsl,'W',W,'maxOLIt',maxOLIt,'P0',POS,'findWindowSize',false);
    tic;
    [KFH,PFH] = kalmanFiniteHorizonLTI(A,C,Q,R,E,optsFH);
    toc;
    if ~isnan(KFH)
        fprintf("Trace FH: %f\n",trace(PFH));
        fprintf("Trace FH/C: %f\n",trace(PFH)/trace(PC));
        fprintf("Maximum absolute CL eigenvalue: %f\n",max(abs(eig((eye(size(A,1))-KFH*C)*A))));
        save(sprintf("./models/model%d/synth_FH.mat",seed),'KFH','PFH','optsFH');
    end
end