function flag = stability_mfh(model,W)

    % Load model 
    model_path = sprintf("./models/model%d/model.mat",model);
    load(model_path,'A','C','Q','R','E')

    % Load MFH synthesis
    load(sprintf("./models/model%d/synth_MFH_%d.mat",model,W),'KMFH');

    % Compute matrix according to condition in Section II.B
    Phi = eye(size(A));
    for i = 1:W
        Phi = (eye(size(A))-KMFH{i,1}*C)*A*Phi;
    end
    
    % Check whether it is Shur
    if max(abs(eig(Phi))) >= 1
        flag = 0;
        fprintf("MFH synthesis is unstable.\n");
    else
        flag = 1;
        fprintf("MFH synthesis is stable.\n");
    end

end