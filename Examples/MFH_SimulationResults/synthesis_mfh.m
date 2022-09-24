function [KMFH, PMFH, PMFH_seq] = synthesis_mfh(seed,epsl_inf,W,maxIt)  

    % Load model 
    model_path = sprintf("./models/model%d/model.mat",seed);
    load(model_path,'A','C','Q','R','E')

    % Load centralized performance
    load(sprintf("./models/model%d/synth_C.mat",seed),'PC');
    
    % MFH
    tic;
    optsMFH = struct('verbose',true,'epsl_inf',epsl_inf,'epsl',epsl_inf/10,'maxIt',maxIt);
    [KMFH,PMFH,PMFH_seq] = MHEMovingFiniteHorizonLTI(A,C,Q,R,E,W,optsMFH);
    comp_time = toc;
    fprintf("Trace MFH: %f\n",trace(PMFH));
    fprintf("Trace MFH/C: %f\n",trace(PMFH)/trace(PC));
    save(sprintf("./models/model%d/synth_MFH_%d.mat",seed,W),'KMFH','PMFH','PMFH_seq','W','optsMFH');
    
    % Log
    fprintf("For W_ss = %d (tr = %g) elapsed time is %g seconds.\n", W, trace(PMFH), comp_time);
    
end