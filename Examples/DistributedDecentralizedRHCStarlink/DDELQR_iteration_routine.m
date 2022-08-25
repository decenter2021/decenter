%% Input:   t  
%           x_t dim: 7N x 1
% Output:   MPC_u dim: 3 x N

%% Check redefinition of anchor 
% Compute mean orbital elements
OE_t = zeros(n_single,N);
%for i = 1:N
parfor i = 1:N
    % With package osculating2mean
    aux = rv2OEOsc(x_t((i-1)*(n_single+1)+1:i*(n_single+1)-1));
    OE_t(:,i) = OEOsc2OEMeanEUK(t,aux,degree_osculating2mean);
end
% Condition for nominal anchor redefinition
% Define nominal anchor just by the initial state
if t < Tctrl/2 % Did not use t < 0 to avoid numerical issues
    fprintf('@MATLAB controller: Redefining nominal anchor at t = %g s.\n',t);
    [u0,Omega0] = nominalConstellationAnchor(OE_t([2 6],:),walkerParameters);
    t0 = t;
    walkerAnchor = [t0; u0; Omega0];
end

%% Transformation to relative obit elements
dalpha_t = zeros(n_single,N);    
parfor i = 1:N
%for i = 1:N
    dalpha_t(:,i) = OEMean2OERel(OE_t(:,i),...
         nominalConstellationOEMean(t,i,walkerParameters,semiMajorAxis,walkerAnchor,0));
end

%% Check if a new MPC window needs to be computed
% Keep using the current window
if MPC_currentWindowGain < MPC_d && t > Tctrl/2
    MPC_currentWindowGain = MPC_currentWindowGain+1;
% Compute a new window and use the first gain
else
    MPC_currentWindowGain = 1;  
    %% Compute new MPC window
    fprintf('@MATLAB controller: Computing new MPC window starting at t = %g s.\n',t);
    
    %% MPC centralized backwards loop 
    % Initializations of dynamics matrices and topology for t = MPC_T
    tau = MPC_T;
    % Time instant (s) of MPC iteration
    t_tau = t+tau*Tctrl;

    % Predict nominal contellation at MPC iteration
    OEMeanNominal_tau = nominalConstellationOEMean(t_tau,1:N,walkerParameters,semiMajorAxis,walkerAnchor,0);  
    % Scheduled time of communication regarding this instant
    t_com = t - (tau+2)*Tt;
    OEMeanNominal_t_com = nominalConstellationOEMean(t_com,1:N,walkerParameters,semiMajorAxis,walkerAnchor,0);
    xNominal_tau = zeros(n_single,N);
    xNominal_t_com = zeros(n_single,N);
    parfor i = 1:N
        xNominal_tau(:,i) = OEOsc2rv(OEMeanEU2OEOsc(OEMeanNominal_tau(:,i)));
        xNominal_t_com(:,i) = OEOsc2rv(OEMeanEU2OEOsc(OEMeanNominal_t_com(:,i)));
    end

    % Predict topology
    % With communication restrictions  
    % Compute in-neighbourhood
    parfor i = 1:N            
        % If there were no restrictions
        aux = LEOConstellationTrackingGraph_DMinus(i,xNominal_tau,trackingRange,trackingmaxInNeighbourhood); 
        % Communication restriction mask
        mask = false(size(aux));
        for l = 1:length(aux)
            if norm(xNominal_t_com(1:3,i)-xNominal_t_com(1:3,aux(l))) < LOSRange
                mask(l) = true;
            end
        end
        % Apply mask
        Di_tau_minus{i,1} = aux(mask);
    end
    
    % Without communication restrictions
%     % Compute in-neighbourhood
%     parfor i = 1:N 
%         Di_tau_minus{i,1} = LEOConstellationTrackingGraph_DMinus(i,xNominal_tau,trackingRange,trackingmaxInNeighbourhood); 
%     end
    
    % Compute out-neighbourhood
    parfor i = 1:N
        Di_tau{i,1} = LEOConstellationTrackingGraph_DPlus(i,xNominal_tau,trackingRange,Di_tau_minus); 
    end

    % Predict output dynamics and state weights
    %for i = 1:N
    parfor i = 1:N
        % Predict output dynamics and state weights
        [Haux,Q{i,1},~,o(i)] = LEOConstellationTrackingDynamics(i,Di_tau_minus,inclination);
        aux = 1;
        % Temporary variable to allow paralelization of cell H
        H_tmp = cell(1,N);
        for j = Di_tau_minus{i,1}'
            H_tmp{1,j} = Haux{aux,1};
            aux = aux + 1;
        end
        H(i,:) = H_tmp;
    end
    
    % MPC iterations start only at t = MPC_T-1
    for tau = MPC_T-1:-1:0
        %% Step 1. Predict topology and distributed dynamics       
        % Time instant (s) of MPC iteration
        t_tau = t+tau*Tctrl;
        %% 1.1 Predict nominal constellation
        % Predict nominal contellation at MPC iteration
        OEMeanNominal_tau = nominalConstellationOEMean(t_tau,1:N,walkerParameters,semiMajorAxis,walkerAnchor,0);   
        % Scheduled time of communication regarding this instant
        t_com = t - (tau+2)*Tt;
        OEMeanNominal_t_com = nominalConstellationOEMean(t_com,1:N,walkerParameters,semiMajorAxis,walkerAnchor,0);
        xNominal_tau = zeros(n_single,N);
        xNominal_t_com = zeros(n_single,N);
       
        parfor i = 1:N
            xNominal_tau(:,i) = OEOsc2rv(OEMeanEU2OEOsc(OEMeanNominal_tau(:,i)));
            xNominal_t_com(:,i) = OEOsc2rv(OEMeanEU2OEOsc(OEMeanNominal_t_com(:,i)));
        end
        %% 1.2 Predict tracking topology for tau
        % Store topology of last iteration (tau+1)
        Di_tau_1 = Di_tau;

        % With communication restrictions
        % Compute in-neighbourhood
        parfor i = 1:N            
            % If there were no restrictions
            aux = LEOConstellationTrackingGraph_DMinus(i,xNominal_tau,trackingRange,trackingmaxInNeighbourhood); 
            % Communication restriction mask
            mask = false(size(aux));
            for l = 1:length(aux)
                if norm(xNominal_t_com(1:3,i)-xNominal_t_com(1:3,aux(l))) < LOSRange
                    mask(l) = true;
                end
            end
            % Apply mask
            Di_tau_minus{i,1} = aux(mask);
        end
        
        % Without communication restrictions
        % Compute in-neighbourhood
%         parfor i = 1:N 
%             Di_tau_minus{i,1} = LEOConstellationTrackingGraph_DMinus(i,xNominal_tau,trackingRange,trackingmaxInNeighbourhood); 
%         end
        
        % Compute out-neighbourhood
        parfor i = 1:N
            Di_tau{i,1} = LEOConstellationTrackingGraph_DPlus(i,xNominal_tau,trackingRange,Di_tau_minus); 
        end
        
        
        %% Step 2: Recompute covarinaces (Compute P(tau+1) )
        if decentralized       
            P_kl_prev = P_kl;
            if tau ~= MPC_T-1
                % Computations done distributedly across agents
                %for i = 1:N
                parfor i = 1:N
                    P_kl{i,1} = newCovarianceStorage(Di_tau{i,1});
                    % Compute each block of P
                    for j = 1:size(P_kl{i,1},1)     
                        p = P_kl{i,1}{j,1}(1);
                        q = P_kl{i,1}{j,1}(2);
                        Cp = zeros(length(Di_tau{p,1})*n_single,n_single);
                        Cq = zeros(length(Di_tau{q,1})*n_single,n_single);
                        Prs = zeros(length(Di_tau{p,1})*n_single,length(Di_tau{q,1})*n_single);
                        lossPrs = 0;
                        count_r = 0; 
                        % Sum over indices r ans s
                        % Build matrix P_rs, C_p, and Cq for each (p,q)
                        for r = Di_tau_1{p,1}'    
                            if p == q
                                count_s = count_r;
                                for s = Di_tau_1{q,1}(count_r+1:end)'
                                    % Get the P_rs computation available to i
                                    [aux,loss] = searchP(i,r,s,P_kl_prev,Di_tau_1);
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
                                for s = Di_tau_1{q,1}'
                                    % Get the P_rs computation available to i
                                    [aux,loss] = searchP(i,r,s,P_kl_prev,Di_tau_1);
                                    Prs(count_r*n_single+1:(count_r+1)*n_single,count_s*n_single+1:(count_s+1)*n_single)=...
                                        aux;             
                                    if count_r == 0
                                        %Cq(count_s*n_single+1:(count_s+1)*n_single,:) = A{q,1}*(q==s)-B{s,1}*getK(MPC_K{s,tau+2},q);   
                                        Cq(count_s*n_single+1:(count_s+1)*n_single,:) = A{q,1}*(q==s)-B{s,1}*K_tau_1{q,1}{count_s+1,1};   
                                    end 
                                    lossPrs = lossPrs + loss;
                                    count_s = count_s + 1;    
                                end
                            end                 
                            Cp(count_r*n_single+1:(count_r+1)*n_single,:) = A{p,1}*(p==r)-B{r,1}*K_tau_1{p,1}{count_r+1,1};      
                            count_r = count_r + 1;
                        end               
                       
                        % If p == q then, P_rs should be positive definite
                        if p == q
                            % If p == q then, P_rs should be positive definite
                            Prs = forcePositiveDefiniteness(Prs);
                            Cq = Cp;
                        end
                        % Intersection of D_p+ and D_q+
                        [r_cap_pq,r_cap_pq_idxp,r_cap_pq_idxq] = intersect(Di_tau_1{p,1},Di_tau_1{q,1});
                        aux = zeros(n_single,n_single);
                        for l = 1:length(r_cap_pq)
                            % parfor does not allow to assign a variable
                            % named 'r'
                            r_ = r_cap_pq(l);
                            r_idxp = r_cap_pq_idxp(l);
                            r_idxq = r_cap_pq_idxq(l);
                            aux = aux + H{r_,p}'*Q{r_,1}*H{r_,q} + K_tau_1{p,1}{r_idxp,1}'*R{r_,1}*K_tau_1{q,1}{r_idxq,1};
                                    % + getK(MPC_K{r_,tau+2},p)'*R{r_,1}*getK(MPC_K{r_,tau+2},q);
                        end
                        P_kl{i,1}{j,2} = Cp'*Prs*Cq + aux;
                        % Update loss
                        P_kl{i,1}{j,3} = lossPrs;
                    end 
                end
            else
                % Computations done distributedly across agents
                %for i = 1:N
                parfor i = 1:N
                    P_kl{i,1} = newCovarianceStorage(Di_tau{i,1});
                    % Compute each block of P
                    for j = 1:size(P_kl{i,1},1)         
                        p = P_kl{i,1}{j,1}(1);
                        q = P_kl{i,1}{j,1}(2);                
                        % Intersection of D_p+ and D_q+
                        [r_cap_pq,~] = intersect(Di_tau_1{p,1},Di_tau_1{q,1});
                        aux = zeros(n_single,n_single);                    
                        for r = r_cap_pq'
                            aux= aux + H{r,p}'*Q{r,1}*H{r,q};
                        end
                        P_kl{i,1}{j,2} = aux;
                        % Update loss (not necessary, it is set to 0 on creation)
                        % P_kl{i,1}{j,3} = 0;
                    end
                end
            end
        else
            % Centralized version below:
            % 2.1.1 Compute global matrices (H(tau+1), Q(tau+1), R(tau+1), K(tau+1), B(tau+1), A(tau+1))
            % Init matrices 
            og = sum(o)*o_single_rel+N*o_single_self;              
            Hg = zeros(og,N*n_single); 
            Qg = zeros(og);
            % Build matrices
            aux = 0;       
            for i = 1:N     
                for j = Di_tau_1{i,1}'
                    Hg(aux+1:aux+o(i)*o_single_rel + o_single_self,(j-1)*n_single+1:(j-1)*n_single+n_single) = H{i,j};  
                end
                Qg(aux+1:aux+o(i)*o_single_rel + o_single_self, aux+1:aux+o(i)*o_single_rel + o_single_self) = Q{i,1};         
                aux = aux+o(i)*o_single_rel + o_single_self;
            end
            % Compute P(tau+1)
            if tau ~= MPC_T-1
                P_tau_1 = Hg'*Qg*Hg + Kg'*Rg*Kg + (Ag-Bg*Kg)'*P_tau_1*(Ag-Bg*Kg);
            else
                P_tau_1 = Hg'*Qg*Hg;
            end
        end
        
        %% Step 3: Predict dynamic and output dynamics matricexs for tau (A(tau), B(tau), R(tau))
        %for i = 1:N
        parfor i = 1:N
            % Predict  dynamics
            [A{i,1},B{i,1}] = STMSatellite(OEMeanNominal_tau(:,i),Tctrl);
            % Predict output dynamics, state weights, and input weights
            [Haux,Q{i,1},R{i,1},o(i)] = LEOConstellationTrackingDynamics(i,Di_tau_minus,inclination);
            aux = 1;
            % Temporary variable to allow paralelization of cell H
            H_tmp = cell(1,N);
            for j = Di_tau_minus{i,1}'
                H_tmp{1,j} = Haux{aux,1};
                aux = aux + 1;
            end
            H(i,:) = H_tmp;
        end
        
        %% Step 4: Compute gains
        if decentralized
            % Decentralized gain computation 
            %for i = 1:N
            parfor i = 1:N
                % Compute augmented innovation covariance matrix
                Si = zeros(m_single*length(Di_tau{i,1}));
                % Compute augmented P_i tilde
                Pi = zeros(m_single*length(Di_tau{i,1}),n_single);
                count = 1;
                for pidx = 1:length(Di_tau{i,1})
                    for qidx = pidx:length(Di_tau{i,1})
                        p = Di_tau{i,1}(pidx);
                        q = Di_tau{i,1}(qidx);                   
                        % Fill Si
                        Si((pidx-1)*m_single+1:(pidx-1)*m_single+m_single,(qidx-1)*m_single+1:(qidx-1)*m_single+m_single) = ...
                            B{p,1}'*P_kl{i,1}{count,2}*B{q,1} + (p==q)*R{p,1};
                        if pidx ~= qidx
                            Si((qidx-1)*m_single+1:(qidx-1)*m_single+m_single,(pidx-1)*m_single+1:(pidx-1)*m_single+m_single) = ...
                                Si((pidx-1)*m_single+1:(pidx-1)*m_single+m_single,(qidx-1)*m_single+1:(qidx-1)*m_single+m_single)';
                        end
                        % Fill Pi
                        if p == i || q == i
                            if p == i
                                Pi((qidx-1)*m_single+1:(qidx-1)*m_single+m_single,:) = B{q,1}'*P_kl{i,1}{count,2}'*A{i,1};
                            else
                                Pi((pidx-1)*m_single+1:(pidx-1)*m_single+m_single,:) = B{p,1}'*P_kl{i,1}{count,2}*A{i,1};
                            end
                        end
                        count = count + 1;
                    end                  
                end
                % Compute gains 
                Ki = Si\Pi;
                % Fill K_tau_1
                K_tau_1{i,1} = cell(length(Di_tau{i,1}),1);
                for pidx = 1:length(Di_tau{i,1})
                    K_tau_1{i,1}{pidx,1} = Ki((pidx-1)*m_single+1:(pidx-1)*m_single+m_single,:);
                end
            end
        else
            % Centralized computation below
            % 4.1 Compute global matrices (B, R, A)
            Bg = zeros(n_single*N,m_single*N);
            Rg = zeros(m_single*N);
            Ag = zeros(n_single*N);      
            % Predict dynamics
            for i = 1:N
                Ag((i-1)*n_single+1:i*n_single,(i-1)*n_single+1:i*n_single) = A{i,1};
                Bg((i-1)*n_single+1:i*n_single,(i-1)*m_single+1:i*m_single) = B{i,1};
                Rg((i-1)*m_single+1:i*m_single,(i-1)*m_single+1:i*m_single) = R{i,1};
            end
            % 4.2 Compute gains
            Sg = Bg'*P_tau_1*Bg + Rg;
            Kg = Sg\Bg'*P_tau_1*Ag;
        end
        
        %% Step 3.5 Transmit gains
        if decentralized
            % Fill decentralized MPC gains k(tau)
            if tau <= MPC_d    
                % If the topology is undirected the implementation below is
                % more efficient 
%                 for i = 1:N
%                 %parfor i = 1:N
%                     MPC_K{i,tau+1} = cell(length(Di_tau{i,1}),2); % MPC_K{i,tau+1} = cell(length(Di_tau_minus{i,1}),2);
%                     % They gains must be recived by agents in D_i-, but given
%                     % that, in this case in particular, the edges are
%                     % undirected, then we can use D_i+ (Di_tau) to index the
%                     % agents and retrieve the gains (it is more efficient)
%                     % Warning: if the topology is directed the code above must
%                     % be adapted (we have to compute and use Di_tau_minus)
%                     for pidx = 1:length(Di_tau{i,1}) % pidx = 1:length(Di_tau_minus{i,1})
%                         p = Di_tau{i,1}(pidx); % p = Di_tau_minus{i,1}(pidx);
%                         MPC_K{i,tau+1}{pidx,1} = p;
%                         i_idx = find(Di_tau{p,1}==i);
%                         MPC_K{i,tau+1}{pidx,2} = K_tau_1{p,1}{i_idx,1};
%                     end
%                 end

                % If the topology is directed, then the implementation
                % above is not valid
                for i = 1:N
                %parfor i = 1:N
                    MPC_K{i,tau+1} = cell(length(Di_tau_minus{i,1}),2);
                    for pidx = 1:length(Di_tau_minus{i,1})
                        p = Di_tau_minus{i,1}(pidx);
                        MPC_K{i,tau+1}{pidx,1} = p;
                        i_idx = find(Di_tau{p,1}==i);
                        MPC_K{i,tau+1}{pidx,2} = K_tau_1{p,1}{i_idx,1};
                    end
                end  
            end
        else       
            % Centralized below
    
%             % For debug: fill K_tau_1 gains for the computation of P_tau_1 in the next MPC
%             % iteration
%             for i = 1:N
%                 K_tau_1{i,1} = cell(length(Di_tau{i,1}),1);
%                 for j = 1:length(Di_tau{i,1})
%                     p = Di_tau{i,1}(j);
%                     K_tau_1{i,1}{j,1} = Kg((p-1)*m_single+1:(p-1)*m_single+m_single,(i-1)*n_single+1:(i-1)*n_single+n_single);
%                 end
%             end
%     
%             % For debug: enforce saprsity of global gain
%             for i = 1:N
%                 for j = 1:N
%                     if sum(i==Di_tau{j,1}) == 0
%                         Kg((i-1)*m_single+1:(i-1)*m_single+m_single,(j-1)*n_single+1:(j-1)*n_single+n_single) = zeros(m_single,n_single);
%                     end
%                 end
%             end
        
            % Fill centralized MPC gains k(tau)
            if tau <= MPC_d           
                for i = 1:N
                    MPC_K{i,tau+1} = cell(N,2);
                    for j = 1:N
                        MPC_K{i,tau+1}{j,1} = j;
                        MPC_K{i,tau+1}{j,2} = Kg((i-1)*m_single+1:i*m_single,(j-1)*n_single+1:j*n_single);  
                    end
                end         
            end
        end    
    end
    fprintf('@MATLAB controller: Computed new MPC window starting at t = %g s.\n',t);
end

%% Compute actuation
% Each satellite computes its actuation based on the known gains
%for i = 1:N
parfor i = 1:N
    % Auxiliary varibale to allow parallelization
    uaux = zeros(m_single,1);
    for j = 1:size(MPC_K{i,MPC_currentWindowGain},1)
        uaux = uaux - MPC_K{i,MPC_currentWindowGain}{j,2}*...
            dalpha_t(:,MPC_K{i,MPC_currentWindowGain}{j,1});
    end
    % Compute actuation force from actuation mass
    uaux = uaux*x_t(i*(n_single+1));
    % Saturate thrust
    uaux(uaux>Ct1) = Ct1;
    uaux(uaux<-Ct1) = -Ct1;
    % Assign thrust
    MPC_u(:,i) = uaux;
end


%% Auxiliary functions

%% Auxiliary functins - topology
% Compute tracking graph
% Warning! The output vector must be sorted
% function Di = LEOConstellationTrackingGraph(i,x_hat,separation) 
%     Di = [];
%     for j = 1:size(x_hat,2)
%         if norm(x_hat(1:3,j)-x_hat(1:3,i)) < separation
%             Di = [Di;j];
%         end
%     end
%     % Di = sort(Di); (not necessary in this case)
% end

% Compute tracking graph in-neighbourhood of agent i 
% Warning! The output vector must be sorted
function DPlusi = LEOConstellationTrackingGraph_DPlus(i,x_hat,separation,DMinus) 
    DPlusi = [];
    for j = 1:size(x_hat,2)
        if norm(x_hat(1:3,j)-x_hat(1:3,i)) < separation
            DPlusi = [DPlusi;j];
        end
    end
    % Check if for each agent j in DPlusi, i is in DMinusj
    mask = false(size(DPlusi));
    for j = 1:length(DPlusi)
        mask(j) = sum(DMinus{DPlusi(j),1}==i)>0;
    end
    DPlusi = DPlusi(mask);
    % Sort
    % DPlusi = sort(DPlusi); (not necessary in this case)
end

% Compute tracking graph in-neighbourhood of agent i 
% Warning! The output vector must be sorted
function DMinusi = LEOConstellationTrackingGraph_DMinus(i,x_hat,separation,maxInNeighbourhood) 
    DMinusi = [];
    for j = 1:size(x_hat,2)
        aux = norm(x_hat(1:3,j)-x_hat(1:3,i));
        if aux < separation
            DMinusi = [DMinusi;[j aux]];
        end
    end
    % Check if the number of agents exceeds the limit 
    if size(DMinusi,1) > maxInNeighbourhood+1
        [~,idx] = sortrows(DMinusi,2);
        DMinusi = DMinusi(idx(1:maxInNeighbourhood+1),1);
        DMinusi = sort(DMinusi);
    else
        DMinusi = DMinusi(:,1);
    end
    
end

%% Auxiliary functions - output dynamics
% Compute single satellite output dynamics
function [H,Q,R,o] = LEOConstellationTrackingDynamics(i,Fim,iNominal)
    % Buiding blocks of dynamics 
    % 'Inertial' reference - a,ex,ey,i track to nominal
    H_single_self = [1 zeros(1,5);...
              zeros(3,2) eye(3) zeros(3,1)];
    %Q_single_self = diag(1./([1e-3, 1e-3, 1e-3, 5e-4].^2));
    Q_single_self = diag(1./([1e-4, 0.5e-2, 0.5e-2, 1e-2].^2));

    o_single_self = 4;
    
    H_single_rel = [0 1 zeros(1,3) -1/tan(iNominal);...
             zeros(1,5) 1/sin(iNominal)];
    Q_single_rel = diag(1./([1e-4 1e-4].^2));
    o_single_rel = 2;
    n_single = 6;
    
    gamma = 1e4;
    %gamma = 1e6;
    R_single = gamma*diag(1./([0.068, 0.068, 0.068].^2));
    
    % Compute number of output signals
    o = size(Fim{i,1},1)-1;
    
    % Init output variable C
    H = cell(size(Fim{i,1},1),1);     
    
    % Find index of i 
    for j = 1:size(Fim{i,1},1)
        if(Fim{i,1}(j)==i)
            j_Fim_i = j;
            % Init C
            H{j,1} = zeros(o*o_single_rel+o_single_self,n_single);
        end
    end

    % Counter for number of the local output
    aux = 1;
    % Fill C{i,j} with j ~= i
    for j = 1:size(Fim{i,1},1)
        if(Fim{i,1}(j)==i)
            H{j,1}(o*o_single_rel+1:o*o_single_rel+o_single_self,:) = H_single_self;
            continue;
        end

        % Init C
        H{j,1} = zeros(o*o_single_rel+o_single_self,n_single);        
        H{j,1}((aux-1)*o_single_rel+1:aux*o_single_rel,:) = -H_single_rel;
        H{j_Fim_i,1}((aux-1)*o_single_rel+1:aux*o_single_rel,:) = H_single_rel;
        aux = aux+1;                
    end

    % Init output variable R
    Q = blkdiag(kron(eye(o),Q_single_rel),Q_single_self);
    R = R_single;
end

%% Auxiliary functions - Covariance storage and management

% Get feeback gain of agent s in relation to agent q
function K_sq = getK(K_s,q)
    for j = 1:size(Ks,1)
        if K_s{j,1} == q
            K_sq = Ks{j,2};
        end
    end
end

% These function were copied (or slightly adapted) from DDEK for
% constellations

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

function [aux,loss] = searchP(i,r,s,P_kl_pred,Fim)
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
        sz = nchoosek(size(Fim,1),2)+size(Fim,1);
    else
        sz = size(Fim,1);
    end

    % All combinations = distinct + repeated
    idx = zeros(sz,2);
    % Combinations of nodes in F_i^-
    aux = combnk(Fim,2);
    for j = 1:sz
        % Repeated combinations
        if j <= size(Fim, 1)
            idx(j,:) = Fim(j)*ones(1,2);
        % Combinations of nodes in F_i^-
        else
            idx(j,:) = aux(j-size(Fim,1),:);
        end
    end
    idx = sortrows(idx);
    % Create covariance cell
    P_kl = cell(sz,3);
    % Repeated combinations
    for j = 1:sz
        P_kl{j,1} = idx(j,:)';
        P_kl{j,3} = 0;
    end
end

function P = forcePositiveDefiniteness(P)
    P = (1/2)*(P+P');
    [V,D] = eig(P);
    d = diag(D);
    dmin = min(d(d>eps));
    d(d<eps) = max([dmin eps*1e6]);
    D = diag(d);
    P = V*D/V;
end
