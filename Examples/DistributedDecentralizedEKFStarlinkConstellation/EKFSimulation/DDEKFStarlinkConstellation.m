%% Simulation of constellation of LEO satellites with Cartesian coodinates
%% Init 
clear; % Clear workspace variables
tic; % Log the execution time

%% Define constellation
numberOfPlanes = 72;
numberOfSatellitesPerPlane = 22;
% Maximum communication distance normalized by the arc length between
% staellites on the same orbit
ISLRange = 750e3;
semiMajorAxis = 6921000;
% Number of satellites

N = numberOfPlanes*numberOfSatellitesPerPlane;

% Init dimensions of the dynamics of each satellite
o_single = 3;
n_single = 6;
fprintf("Constellation defined.\n");

%% Simulation
Ts = 1; % Sampling time (s)
Tsim = 10;%5730;
ItSim = Tsim/Ts+1; % One week
% Load true data 
load('./data/output_2022_01_22.mat','x');

% Reduce dimension of imported data arrays
parfor i = 1:N
    x{i,1} = x{i,1}(:,1:ItSim);
end
fprintf("Constellation simulation uploaded.\n");

%% Filter simulation - variable definition
%%%% DEK simulation variables
% Estimate vector time series
x_hat = cell(N,1);
for i = 1:N
    x_hat{i,1} = zeros(n_single,ItSim);
end
x_hat_pred = cell(N,1);
% Output vector
y = cell(N,1);
% In neighbourhood
Fim = cell(N,1);
FimNew = cell(N,1);
FimHist = cell(ItSim-1,1);
% Access to inertial measurements (skip = 1)
Fii = [];
% count = 0;
% for i = 0:numberOfPlanes-1
%     j = i*numberOfSatellitesPerPlane+1+count;
%     count = count + 1;
%     count = rem(count,numberOfSatellitesPerPlane);
%     Fii = [Fii j];
% end
Fii = (1:N)';

% Data structures of the dynamics to emulate communication
A = cell(N,1);
Q = cell(N,1);
C = cell(N,N);
o = zeros(N,1);
R = cell(N,N);
% Covariance between nodes in F_i^-
P_kl = cell(N,1);
P_kl_pred = cell(N,1);
% Matrices S_ii
S = cell(N,1);
% Gains K_i
K = cell(N,1);

%% Filter simulation - covariance initializaton
%%%% Covariance initialization
% Initial estimation error covariance
P0_single = blkdiag(10^2*eye(3),0.1^2*eye(3));
%%%% Estimate initialization
for i = 1:N
    x_hat{i,1}(:,1) = x{i,1}(1:n_single,1)+ mvnrnd(zeros(n_single,1),P0_single)';
end

%% Evolution of feedback variables
trace_log = zeros(N,ItSim-1);
P_pos_log = zeros(N,ItSim-1,3);

%% Filter simulation - filter iterations
fprintf("Simulating DEKF.\n");
fprintf("Iteration: %08d/%08d",1,ItSim);

for t = 1:ItSim-1
    % Update counter
    fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%08d/%08d",t,ItSim-1);
    % To distribute less memory across cores define temporary variables
    parfor i = 1:N
        x_t1(:,i) = x{i,1}(1:n_single,t+1);
    end
    %% DEKF - Steps for each time-step for each satellite
    FimHist{t,1} = Fim;
    %% Step 1 - Predict the estimate, topology and compute linearized dynamics 
%     % Before UKF
%     parfor i = 1:N
%         % x_hat_i(t+1|t)
%         % A_ii(t)
%         % Q_ii(t)
%         [x_hat_pred{i,1},A{i,1},Q{i,1}] = predictionCart(x_hat{i,1}(:,t),zeros(3,1),Ts);
%     end
    parfor i = 1:N
        Q{i,1} = QCartUT(x_hat{i,1}(:,t),zeros(3,1),Ts);        
    end
      
    %% Step 2 - Predict the topology 
    % Here we are using the real distance between satellites 
    % We could, instead use an estimate, and if the satellite does not
    % communicate then we disregard it in function
    % LEOConstellationMeasurementGraphSynch
    aux = zeros(n_single,N);
    parfor i = 1:N
        aux(:,i) = x_t1(:,i);
    end
    parfor i = 1:N
        FimNew{i,1} = LEOConstellationMeasurementGraph(i,aux,ISLRange); 
    end

    %% Step 3,4 - Communication
    %% Step 5 - Topology synch
    FimNew = LEOConstellationMeasurementGraphSynch(FimNew);
    %% Step 6 - Update filtered covariance
    % P(t|t)
    if t > 1
        parfor i = 1:N
            P_kl{i,1} = newCovarianceStorage(FimNew{i,1});
            for j = 1:size(P_kl{i,1},1)         
                p = P_kl{i,1}{j,1}(1);
                q = P_kl{i,1}{j,1}(2);
                Cp = zeros(n_single,length(FimNew{p,1})*n_single);
                Cq = zeros(n_single,length(FimNew{q,1})*n_single);
                Prs = zeros(length(FimNew{p,1})*n_single,length(FimNew{q,1})*n_single);
                lossPrs = 0;
                count_r = 0;
                % Sum over indices r ans s
                % Build matrix P_rs, C_p, and Cq for each (p,q)
                for r = Fim{p,1}'    
                    if p == q
                        count_s = count_r;
                        for s = Fim{q,1}(count_r+1:end)'
                            % Get the P_rs computation available to i
                            [aux,loss] = searchP(i,p,q,r,s,P_kl_pred,Fim);
                            Prs(count_r*n_single+1:(count_r+1)*n_single,count_s*n_single+1:(count_s+1)*n_single)=...
                                aux;  
                            if r ~= s
                                Prs(count_s*n_single+1:(count_s+1)*n_single,count_r*n_single+1:(count_r+1)*n_single)=...
                                    Prs(count_r*n_single+1:(count_r+1)*n_single,count_s*n_single+1:(count_s+1)*n_single)';
                            end
                            lossPrs = lossPrs + loss;
                            count_s = count_s + 1;    
                        end   
                    else
                        count_s = 0;
                        for s = Fim{q,1}'
                            % Get the P_rs computation available to i
                            [aux,loss] = searchP(i,p,q,r,s,P_kl_pred,Fim);
                            Prs(count_r*n_single+1:(count_r+1)*n_single,count_s*n_single+1:(count_s+1)*n_single)=...
                                aux;             
                            if count_r == 0
                                Cq(:,count_s*n_single+1:(count_s+1)*n_single) = eye(n_single)*(q==s)-K{q,1}*C{q,s};     
                            end 
                            lossPrs = lossPrs + loss;
                            count_s = count_s + 1;    
                        end
                    end                 
                    Cp(:,count_r*n_single+1:(count_r+1)*n_single) = eye(n_single)*(p==r)-K{p,1}*C{p,r};     
                    count_r = count_r + 1;
                end
               
                if p == q
                    % If p == q then, P_rs should be positive definite
                    Prs = forcePositiveDefiniteness(Prs,0);
                    P_kl{i,1}{j,2} = K{p,1}*R{p,p}*K{p,1}'+Cp*Prs*Cp';
                    P_kl{i,1}{j,2} = (1/2)*(P_kl{i,1}{j,2}+P_kl{i,1}{j,2}');                   
                elseif sum(q == Fim{p,1}) % if p and q have a relative mesurment between them
                    P_kl{i,1}{j,2} = K{p,1}*R{p,q}*K{q,1}' + Cp*Prs*Cq';
                else
                    P_kl{i,1}{j,2} = Cp*Prs*Cq';
                end  
                % Update loss
                P_kl{i,1}{j,3} = lossPrs;
            end 
        end

    else
        % Covariance initialization
        for i = 1:N
            P_kl{i,1} = newCovarianceStorage(FimNew{i,1});
            % Repeated combinations
            for j = 1:size(P_kl{i,1},1)
                p = P_kl{i,1}{j,1}(1);
                q = P_kl{i,1}{j,1}(2);
                if p == q
                    P_kl{i,1}{j,2} = P0_single;
                else
                    P_kl{i,1}{j,2} = zeros(n_single,n_single); 
                end  
            end
        end
    end
    
    %% Step 7 - Update C and R
    % C_ij(t+1)
    % R_ij(t+1)
    parfor i = 1:N
        %[Caux,R{i,1},o(i)] = LEOConstellationOutputDynamics(i,FimNew{i,1},sum(Fii==i));
        [Caux,Raux,o(i)] = LEOConstellationOutputDynamics(i,FimNew,Fii);
        aux = 1;
        % Temporary variable to allow paralelization of cell C
        C_tmp = cell(1,N);
        for j = FimNew{i,1}'
            C_tmp{1,j} = Caux{aux,1};
            aux = aux + 1;
        end
        C(i,:) = C_tmp;
        % Temporary variable to allow paralelization of cell R
        aux = 1;
        for j = FimNew{i,1}'
            C_tmp{1,j} = Raux{aux,1};
            aux = aux + 1;
        end
        R(i,:) = C_tmp;
    end
      
%     %% Step 8 - Predict covariance (before UKF)
%     % P_i(pq)(t+1|t)
%     % S_ii(t+1)
%     % Ki(t+1)
%     P_kl_pred = P_kl;
%     parfor i = 1:N
%         for j = 1:size(P_kl{i,1},1)
%             p = P_kl_pred{i,1}{j,1}(1);
%             q = P_kl_pred{i,1}{j,1}(2);
%             % Local covariance prediction
%             P_kl_pred{i,1}{j,2} = A{p,1}*P_kl{i,1}{j,2}*A{q,1}';
%             if p == q
%                 % if p == q P_i(pq)(t+1|t) must be positive definite
%                 P_kl_pred{i,1}{j,2} = P_kl_pred{i,1}{j,2} + Q{p,1};
%                 P_kl_pred{i,1}{j,2} = (1/2)*(P_kl_pred{i,1}{j,2}+P_kl_pred{i,1}{j,2}');
%             end
%         end
%         
%         % Compute innovation covarinace
%         Pi = zeros(length(FimNew{i,1})*n_single);
%         Ci = zeros(o_single*o(i),length(FimNew{i,1})*n_single);
%         idxi = 0;
%         countk = 0;
%         for k = FimNew{i,1}'
%             if k == i
%                 idxi = countk;
%             end
%             Ci(:,countk*n_single+1:(countk+1)*n_single) = C{i,k};
%             countl = countk;
%             for l = FimNew{i,1}(countk+1:end)'
%                 Pi(countk*n_single+1:(countk+1)*n_single,countl*n_single+1:(countl+1)*n_single)=...
%                     getP(P_kl_pred{i,1},k,l);
%                 if l ~= k
%                     Pi(countl*n_single+1:(countl+1)*n_single,countk*n_single+1:(countk+1)*n_single)=...
%                         Pi(countk*n_single+1:(countk+1)*n_single,countl*n_single+1:(countl+1)*n_single)';
%                 end
%                 countl = countl + 1;
%             end
%             countk = countk + 1;         
%         end
%         % Pi must be positive definite
%         Pi = forcePositiveDefiniteness(Pi,0);
%         % Local computation of innovation covariance matrix
%         S{i,1} = R{i,i}+Ci*Pi*Ci';  
%         S{i,1} = (1/2)*(S{i,1}+S{i,1}');
%         aux = Pi(idxi*n_single+1:(idxi+1)*n_single,:)*Ci';
%          
%         % Compute the gain locally
%         K{i,1} = aux/S{i,1};
%         
%     end

    %% Step 8 - Predict covariance (UKF)
    % P_i(pq)(t+1|t)
    % S_ii(t+1)
    % Ki(t+1)
    P_kl_pred = P_kl;
    parfor i = 1:N
        % Define local concatenated covariance matrices and local concatenated estimates
        % Take avatage of the loop to compute the concatenation of C
        x_hat_i = zeros(length(FimNew{i,1})*n_single,1);
        Pi = zeros(length(FimNew{i,1})*n_single);
        Ci = zeros(o_single*o(i),length(FimNew{i,1})*n_single);
        idxi = 0;
        countk = 0;
        for k = FimNew{i,1}'
            if k == i
                idxi = countk;
            end
            x_hat_i(countk*n_single+1:(countk+1)*n_single,1) = x_hat{i,1}(:,t);
            Ci(:,countk*n_single+1:(countk+1)*n_single) = C{i,k};
            countl = countk;
            for l = FimNew{i,1}(countk+1:end)'
                Pi(countk*n_single+1:(countk+1)*n_single,countl*n_single+1:(countl+1)*n_single)=...
                    getP(P_kl{i,1},k,l);
                if l ~= k
                    Pi(countl*n_single+1:(countl+1)*n_single,countk*n_single+1:(countk+1)*n_single)=...
                        Pi(countk*n_single+1:(countk+1)*n_single,countl*n_single+1:(countl+1)*n_single)';
                end
                countl = countl + 1;
            end
            countk = countk + 1;         
        end    
        % Pi must be positive definite
        Pi = forcePositiveDefiniteness(Pi,0);
        % Compute sigma points and weights
        % Wan, E.A. and Van Der Merwe, R., 2000, October. The unscented Kalman 
        % filter for nonlinear estimation. In Proceedings of the IEEE 2000 Adaptive
        % Systems for Signal Processing, Communications, and Control Symposium 
        % (Cat. No. 00EX373) (pp. 153-158). Ieee.
        % Julier, S.J. and Uhlmann, J.K., 2004. Unscented filtering and 
        % nonlinear estimation. Proceedings of the IEEE, 92(3), pp.401-422.
        % Unscented transformation parameters
        L = length(x_hat_i);
        UTpar_alpha = 1e-3;
        UTpar_k = 0;
        UTpar_beta = 2;
        % Sigma points
        UTpar_lambda = (UTpar_alpha^2)*(L+UTpar_k)-L;
        P_sqrt = chol((L+UTpar_lambda)*Pi)';    
        Ksi = [x_hat_i x_hat_i+P_sqrt x_hat_i-P_sqrt];
        % Weights
        W_m = [UTpar_lambda/(L+UTpar_lambda); ones(2*L,1)*1/(2*(L+UTpar_lambda))];
        W_c = W_m;
        W_c(1) = W_m(1) + 1-UTpar_alpha^2+UTpar_beta; 
        % Propagate sigma points
        for l = 1:2*L+1
            for s = 1:length(FimNew{i,1})
                Ksi((s-1)*n_single+1:s*n_single,l) = predictionCartUT(Ksi((s-1)*n_single+1:s*n_single,l),zeros(3,1),Ts);
            end
        end
        % Compute mean
        x_hat_i = zeros(L,1);
        for l = 1:2*L+1
            x_hat_i = x_hat_i + W_m(l)*Ksi(:,l);
        end
        % Compute covariance
        Ksi = Ksi - x_hat_i;
        Pi = zeros(L,L);
        for l = 1:2*L+1
            Pi = Pi + W_c(l)*(Ksi(:,l)*Ksi(:,l)');
        end      
        % Fill predicted estimate of the satellite
        x_hat_pred{i,1} = x_hat_i(idxi*n_single+1:(idxi+1)*n_single,1);
        % Add additive process noise and fill data structure with covariances
        countk = 0;
        for k = FimNew{i,1}'  
            Pi(countk*n_single+1:(countk+1)*n_single,countk*n_single+1:(countk+1)*n_single) = ...
                Pi(countk*n_single+1:(countk+1)*n_single,countk*n_single+1:(countk+1)*n_single) + Q{k,1};
            countl = countk;
            for l = FimNew{i,1}(countk+1:end)'
                % Find the position on the data structure
                for s = 1:size(P_kl_pred{i,1})
                    if(sum(P_kl_pred{i,1}{s,1} == [k;l]) == 2)
                        P_kl_pred{i,1}{s,2} = ...
                            Pi(countk*n_single+1:(countk+1)*n_single,countl*n_single+1:(countl+1)*n_single);
                        break;
                    end
                end              
                countl = countl + 1;
            end
            countk = countk + 1;         
        end          
        % Local computation of innovation covariance matrix
        S{i,1} = R{i,i}+Ci*Pi*Ci';  
        S{i,1} = (1/2)*(S{i,1}+S{i,1}');
        aux = Pi(idxi*n_single+1:(idxi+1)*n_single,:)*Ci';
        % Compute the gain locally
        K{i,1} = aux/S{i,1};      
    end
    
    %% Step 7 - Take the measurement 
    % Correlated measurement noise generated locally
    for i = 1:N
        % Add error to each measurement relative to satellites in Fim of i
        y{i,1} = zeros(o(i)*o_single,1);
        % Counter for the number of the local output
        aux = 1;         
        for j = 1:size(FimNew{i,1},1)
            % temporary variable to allow paralelization
            y_tmp = zeros(o(i)*o_single,1);
            % Number of satellite 
            satj = FimNew{i,1}(j);
            if satj == i
                continue;
            end
            if satj > i
                y{i,1}((aux-1)*o_single+1:aux*o_single,1) = mvnrnd(zeros(o_single,1),R{i,i}((aux-1)*o_single+1:aux*o_single,(aux-1)*o_single+1:aux*o_single))';
            else
                % Find the number of the local output in j
                auxj = 1;
                for l = 1:size(FimNew{satj,1},1)
                    if FimNew{satj,1}(l) == i
                        break;
                    elseif FimNew{satj,1}(l) ~= satj
                        auxj = auxj + 1;
                    end             
                end
                y{i,1}((aux-1)*o_single+1:aux*o_single,1) = -y{satj,1}((auxj-1)*o_single+1:auxj*o_single);
            end
            aux = aux + 1;           
        end
        % Add inertial error
        if o(i)==aux % if has inertial measurment
            y{i,1}((aux-1)*o_single+1:aux*o_single,1) = mvnrnd(zeros(o_single,1),R{i,i}((aux-1)*o_single+1:aux*o_single,(aux-1)*o_single+1:aux*o_single))';
        end
    end  
    % Take measurments
    parfor i = 1:N
        for j = FimNew{i,1}'
            y{i,1} = y{i,1} + C{i,j}*x_t1(:,j);
        end
    end
    
    %% Step 8 - Update the estimate 
    parfor i = 1:N
        % Compute predicted output (it is linear in cartesian coodinates)
        y_hat = zeros(size(C{i,i},1),1);
        for j = FimNew{i,1}'
            y_hat = y_hat + C{i,j}*x_hat_pred{j,1};
        end
        % Update step
        x_hat{i,1}(:,t+1) = x_hat_pred{i,1}+K{i,1}*(y{i,1}-y_hat);
        Fim{i,1} = FimNew{i,1};
    end
    
    %% Log estimated trace
    parfor i = 1:N
        aux = getP(P_kl{i,1},i,i);
        trace_log(i,t) = trace(aux);
        for m = 1:3
            P_pos_log(i,t,m) = aux(m,m);
        end
    end
end        
%% Evaluate simulation error 
% Compute estimation error in cartesian estimates
error = cell(N,1);
for i = 1:N
    error{i,1} = zeros(n_single,ItSim);
    for t = 1:ItSim
        error{i,1}(:,t) = x_hat{i,1}(:,t)-x{i,1}(1:n_single,t);
    end 
end
fprintf(".\nSimulation completed. Saving data...\n");
save('xhat_simulation.mat','x_hat','Fii','FimHist','trace_log','error','P_pos_log');   
toc
fprintf("Execution completed with success.\n");

%% Auxiliary functions - covariance storage and management

% Get P_ij from the variables of each node
function P = getP(P,k,l)
    for i = 1:size(P,1)
        if(sum(P{i,1} == [k;l]) == 2)
            P = P{i,2};
            return;
        elseif (sum(P{i,1} == [l;k]) == 2)
            P = P{i,2}';
            return;
        end
    end
end

function [P,loss] = getPwloss(P,k,l)
    for i = 1:size(P,1)
        if(sum(P{i,1} == [k;l]) == 2)
            loss = P{i,3};
            P = P{i,2};           
            return;
        elseif (sum(P{i,1} == [l;k]) == 2)
            loss = P{i,3};
            P = P{i,2}';           
            return;
        end
    end
    %fprintf("Fatal error: P_kl not found.\n");
    %P = nan;    
end

function [aux,loss] = searchP(i,p,q,r,s,P_kl_pred,Fim)
    aux = zeros(size(P_kl_pred{1,1}{1,2},1));
    loss = 0;
    access_list = [];
    for l = Fim{i,1}'
        if sum(Fim{l,1} == r) && sum(Fim{l,1} == s )
            access_list = [access_list l];
        end
    end
    
    if isempty(access_list)
        loss = 1;
        return;
    end
    
    % Get all matrices available and respective losses
    P = cell(length(access_list),1);
    lossP = zeros(1,length(access_list));
    for l = 1:length(access_list)
        [P{l,1}, lossP(l)] = getPwloss(P_kl_pred{access_list(l),1},r,s);
    end
    % Find matrices with minimum losses
    access_list = find(lossP == min(lossP));
    for l = access_list
        aux = aux + P{l,1};
    end
    % Output average of matrices with minimum losses
    aux = aux/length(access_list);
    
end

function P_kl = newCovarianceStorage(Fim)
    % Init memory for each node
    % Compute distinct combinations
    if size(Fim,1) > 1
        aux = nchoosek(size(Fim,1),2);
    else
        aux = 0;
    end
    % All combinations = distinct + repeated
    P_kl = cell(aux+size(Fim,1),3);
    % Repeated combinations
    for j = 1:size(Fim,1)
        P_kl{j,1} = Fim(j)*ones(2,1);
        P_kl{j,3} = 0;
    end
    % Combinations of nodes in F_i^-
    aux = combnk(Fim,2);
    for j = size(Fim,1)+1:size(Fim,1)+size(aux,1)
        P_kl{j,1} = aux(j-size(Fim,1),:)';
        P_kl{j,3} = 0;
    end
end

function P = forcePositiveDefiniteness(P,i)
    P = (1/2)*(P+P');
    [V,D] = eig(P);
    d = diag(D);
    dmin = min(d(d>eps));
    d(d<eps) = max([dmin eps*1e6]);
    D = diag(d);
    P = V*D/V;
end

%% Auxiliary functions - output dynamics
% Compute single satellite output dynamics
function [C,R,o] = LEOConstellationOutputDynamics(i,Fim,Fii)
    % Buiding blocks of dynamics 
    C_single = [eye(3) zeros(3,3)];
    o_single = 3;
    n_single = 6;
    R_single_rel = 0.1^2*eye(3);
    R_single_in = (10.0^2)*eye(3);
    
    % Compute number of output signals
    o = size(Fim{i,1},1)-1 + 1*sum(Fii==i);
     
    % Init output variable C
    C = cell(size(Fim{i,1},1),1); 
    % Counter for number of the local output
    aux = 1;
    for j = 1:size(Fim{i,1},1)
        C{j,1} = zeros(o*o_single,n_single);
        if(Fim{i,1}(j)==i)
            for l = 1:o
                C{j,1}((l-1)*o_single+1:l*o_single,:) = C_single;
            end
        else       
            C{j,1}((aux-1)*o_single+1:aux*o_single,:) = -C_single;
            aux = aux+1;         
        end
    end
    
    % Init output variable R
    R = cell(size(Fim{i,1},1),1); 
    % Counter for number of the local output
    aux = 1;
    for j = 1:size(Fim{i,1},1)
        satj = Fim{i,1}(j);
        oj = size(Fim{satj,1},1)-1 + 1*sum(Fii==satj);
        R{j,1} = zeros(o*o_single,oj*o_single);
        if(Fim{i,1}(j)==i)
            if sum(Fii==i)
                R{j,1} = blkdiag(kron(eye(o-1),R_single_rel),R_single_in);
            else
                R{j,1} = kron(eye(o),R_single_rel);
            end
        else
            % Find the number of the local output in j
            auxj = 1;
            for l = 1:size(Fim{satj,1},1)
                if Fim{satj,1}(l) == i
                    break;
                elseif Fim{satj,1}(l) ~= satj
                    auxj = auxj + 1;
                end             
            end
            R{j,1}((aux-1)*o_single+1:aux*o_single,(auxj-1)*o_single+1:auxj*o_single) = -R_single_rel;
            aux = aux+1;
        end
    end
end

%% Auxiliary functions - topology

% Compute measurment graph
function Fim = LEOConstellationMeasurementGraph(i,x_hat,separation) 
    % In neighborhood 
    Fim = [];
    for j = 1:size(x_hat,2)
        if norm(x_hat(1:3,j)-x_hat(1:3,i)) < separation
            Fim = [Fim;j];
        end
    end
end

function Fim = LEOConstellationMeasurementGraphSynch(Fim)

end


%% Auxiliary functions - Prediction dynamics

function [x_pred] = predictionCartUT(x,u,Ts)
    %% State prediction with J2 perturbation by soving ODE
    t = [0 Ts]; % Set time instants to evaluate x
    options = odeset('RelTol',1e-10); % ODE solver options
    [~,x] = ode45(@(t,x) f(t,x,u),t,x,options); % Solve ODE
    x_pred = x(end,:)'; % Output solution at next time instant   
end

function [Q] = QCartUT(x,u,Ts)
       %% Q
    % Define Hill frame
%     or = x(1:3)/norm(x(1:3));
%     on = cross(x(1:3),x(4:6));
%     on = on/norm(on);
%     ot = cross(on,or);
%     
%     R_Hill2Cart = [or ot on];
%     
%     R = blkdiag(R_Hill2Cart,R_Hill2Cart);
%     
%     QHill = diag(...);
%     Q =  R*QHill*R';

    Q = zeros(6,6);
    Q(1,1) = 1.9667e-3;
    Q(2,2) = 1.9667e-3;
    Q(3,3) = 1.4577e-3;
    Q(4,4) = 3.3821e-5;
    Q(5,5) = 3.3821e-5;
    Q(6,6) = 4.0904e-5;
    Q(1,4) = 2.5152e-4;
    Q(4,1) = Q(1,4);
    Q(2,5) = 2.5152e-4;
    Q(5,2) = Q(2,5);
    Q(3,6) = 2.4239e-4;
    Q(6,3) = Q(3,6);
    Q = 1000*Q;
end

function [x_pred,Phi,Q] = predictionCart(x,u,Ts)
    %% State prediction with J2 perturbation by soving ODE
    t = [0 Ts/2 Ts]; % Set time instants to evaluate x
    options = odeset('RelTol',1e-10); % ODE solver options
    [~,x] = ode45(@(t,x) f(t,x,u),t,x,options); % Solve ODE
    x_pred = x(end,:)'; % Output solution at next time instant
    %% State transition matrix by RK4 integration
    k1 = dPhidt(x(1,:)',u,eye(6));
    k2 = dPhidt(x(2,:)',u,eye(6)+Ts*k1/2);
    k3 = dPhidt(x(2,:)',u,eye(6)+Ts*k2/2);
    k4 = dPhidt(x(3,:)',u,eye(6)+Ts*k3);
    Phi = eye(6) + (1/6)*Ts*(k1+2*k2+2*k3+k4);
    %% Q
    Q = QCartUT(x,u,Ts);        
end

function dxdt = f(t,x,u)
    %% Constants 
    mu = 3.986004418e14; %(m^3 s^-2)
    RE = 6371e3; %(m)
    J2 = 1082.6267e-6;
    %Ct1 = 0.068; % (N)
    Isp = 1640; % (s)
    g0 = 9.81; % (ms^-2)
    fs = 1;
    %% State vector 
    rx = x(1);
    ry = x(2);
    rz = x(3);
    vx = x(4);
    vy = x(5);
    vz = x(6);
    %m = x(7);
    m = 260;
    fdx = u(1);
    fdy = u(2);
    fdz = u(3);
    %% Compute spacecraft dynamics (cf. Mathematica)
    %dxdt = zeros(7,1);
    dxdt = zeros(6,1);
    % Kinematics
    dxdt(1:3) = x(4:6);
    % Dynamics
    % Central body acceleration
    accKep = zeros(3,1);
    accKep(1) = (-1)*mu*rx*(abs(rx)^2+abs(ry)^2+abs(rz)^2)^(-3/2);
    accKep(2) = (-1)*mu*ry*(abs(rx)^2+abs(ry)^2+abs(rz)^2)^(-3/2);
    accKep(3) = (-1)*mu*rz*(abs(rx)^2+abs(ry)^2+abs(rz)^2)^(-3/2);
    % J2 Acceleration
    accJ2 = zeros(3,1);
    accJ2(1) = (15/2)*J2*mu*RE^2*rx*rz^2*(rx^2+ry^2+rz^2)^(-7/2)+(-3/2)*J2*mu* ...
                RE^2*rx*(rx^2+ry^2+rz^2)^(-5/2);
    accJ2(2) = (15/2)*J2*mu*RE^2*ry*rz^2*(rx^2+ry^2+rz^2)^(-7/2)+(-3/2)*J2*mu* ...
                RE^2*ry*(rx^2+ry^2+rz^2)^(-5/2);
    accJ2(3) = (15/2)*J2*mu*RE^2*rz^3*(rx^2+ry^2+rz^2)^(-7/2)+(-9/2)*J2*mu*RE^2* ...
                rz*(rx^2+ry^2+rz^2)^(-5/2);
    % Thruster 
    accT = zeros(3,1);
    accT(1) = fdx*fs/m;
    accT(2) = fdy*fs/m;
    accT(3) = fdz*fs/m;
    % Complete dynamics
    dxdt(4:6) = accKep + accJ2 + accT;
    % Mass dynamics
    %dxdt(7) = -(sum(abs(u)))/(Isp*g0);
end

function Phi_dot = dPhidt(x,u,Phi)
    %% Constants 
    mu = 3.986004418e14; %(m^3 s^-2)
    RE = 6371e3; %(m)
    J2 = 1082.6267e-6;
    Ct1 = 0.068; % (N)
    Isp = 1640; % (s)
    g0 = 9.81; % (ms^-2)
    fs = 1;
    %% State vector 
    rx = x(1);
    ry = x(2);
    rz = x(3);
    vx = x(4);
    vy = x(5);
    vz = x(6);
    %m = x(7);
    m = 260;
    fdx = u(1); % (u in N)
    fdy = u(2);
    fdz = u(3);
    %% Compute F 
    %F = zeros(7);
    F = zeros(6);
    F(1:3,4:6) = eye(3);
    F(4,1) = (-105/2)*J2*mu*RE^2*rx^2*rz^2*(rx^2+ry^2+rz^2)^(-9/2)+(15/2)*J2* ...
            mu*RE^2*rx^2*(rx^2+ry^2+rz^2)^(-7/2)+(15/2)*J2*mu*RE^2*rz^2*(rx^2+ ...
            ry^2+rz^2)^(-7/2)+(-3/2)*J2*mu*RE^2*(rx^2+ry^2+rz^2)^(-5/2)+3*mu* ...
            rx^2*(rx^2+ry^2+rz^2)^(-5/2)+(-1)*mu*(rx^2+ry^2+rz^2)^(-3/2);
    F(4,2) = (-105/2)*J2*mu*RE^2*rx*ry*rz^2*(rx^2+ry^2+rz^2)^(-9/2)+(15/2)*J2* ...
            mu*RE^2*rx*ry*(rx^2+ry^2+rz^2)^(-7/2)+3*mu*rx*ry*(rx^2+ry^2+rz^2) ...
            ^(-5/2);
    F(4,3) = (-105/2)*J2*mu*RE^2*rx*rz^3*(rx^2+ry^2+rz^2)^(-9/2)+(45/2)*J2*mu* ...
            RE^2*rx*rz*(rx^2+ry^2+rz^2)^(-7/2)+3*mu*rx*rz*(rx^2+ry^2+rz^2)^( ...
            -5/2);
    F(5,1) = (-105/2)*J2*mu*RE^2*rx*ry*rz^2*(rx^2+ry^2+rz^2)^(-9/2)+(15/2)*J2* ...
            mu*RE^2*rx*ry*(rx^2+ry^2+rz^2)^(-7/2)+3*mu*rx*ry*(rx^2+ry^2+rz^2) ...
            ^(-5/2);
    F(5,2) = (-105/2)*J2*mu*RE^2*ry^2*rz^2*(rx^2+ry^2+rz^2)^(-9/2)+(15/2)*J2* ...
            mu*RE^2*ry^2*(rx^2+ry^2+rz^2)^(-7/2)+(15/2)*J2*mu*RE^2*rz^2*(rx^2+ ...
            ry^2+rz^2)^(-7/2)+(-3/2)*J2*mu*RE^2*(rx^2+ry^2+rz^2)^(-5/2)+3*mu* ...
            ry^2*(rx^2+ry^2+rz^2)^(-5/2)+(-1)*mu*(rx^2+ry^2+rz^2)^(-3/2);
    F(5,3) = (-105/2)*J2*mu*RE^2*ry*rz^3*(rx^2+ry^2+rz^2)^(-9/2)+(45/2)*J2*mu* ...
            RE^2*ry*rz*(rx^2+ry^2+rz^2)^(-7/2)+3*mu*ry*rz*(rx^2+ry^2+rz^2)^( ...
            -5/2);
    F(6,1) = (-105/2)*J2*mu*RE^2*rx*rz^3*(rx^2+ry^2+rz^2)^(-9/2)+(45/2)*J2*mu* ...
            RE^2*rx*rz*(rx^2+ry^2+rz^2)^(-7/2)+3*mu*rx*rz*(rx^2+ry^2+rz^2)^( ...
            -5/2);
    F(6,2) = (-105/2)*J2*mu*RE^2*ry*rz^3*(rx^2+ry^2+rz^2)^(-9/2)+(45/2)*J2*mu* ...
            RE^2*ry*rz*(rx^2+ry^2+rz^2)^(-7/2)+3*mu*ry*rz*(rx^2+ry^2+rz^2)^( ...
            -5/2);   
    F(6,3) = (-105/2)*J2*mu*RE^2*rz^4*(rx^2+ry^2+rz^2)^(-9/2)+45*J2*mu*RE^2* ...
            rz^2*(rx^2+ry^2+rz^2)^(-7/2)+(-9/2)*J2*mu*RE^2*(rx^2+ry^2+rz^2)^( ...
            -5/2)+3*mu*rz^2*(rx^2+ry^2+rz^2)^(-5/2)+(-1)*mu*(rx^2+ry^2+rz^2)^( ...
            -3/2);
    %F(4:6,7) = (-1)*u*fs*m^(-2);
    Phi_dot = F*Phi;
end
