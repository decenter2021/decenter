%% Description - Traffic Network Control
% This script simulates the application of 4 decentralized control methods
% to the signal control problem of congested urban traffic networks:
% 1. DTUC with decentralized configuration psi
% 2. DTUC with decentralized configuration phi
% 3. D2TUC with decentralized configuration psi
% 4. D2TUC with decentralized configuration phi
% It also simulates two centralized methods as a baseline:
% 1. TUC, as detailled in [2]
% 2. TUC - QPC (centralized version of D2TUC) inspired in [2] 
% The implementation is detailled in [1].

%% Initilalize workspace
clear;

%% Import LTI model of Chania network
model = modelSynthesis("./data/");

%% Modal decomposition as detailed in [1]
% Separation of controllable and uncontrolable modes
ctrbM = ctrb(model.A,model.Bg); % controlability matrix
r = rank(ctrbM);
Z = size(model.A,1);
H = orth(ctrbM);
V = null(H');
W = [H V];
A_hat = W\model.A*W;
Bg_hat = W\model.Bg; 
Bg1_hat = [eye(r) zeros(r,Z-r)]*Bg_hat;
A1_hat = [eye(r) zeros(r,Z-r)]*A_hat*[eye(r);zeros(Z-r,r)];

%% Generate historic actuation as detailed in [1]
gN_DTUC = -(1/model.C)*(Bg1_hat'*Bg1_hat)\Bg1_hat'*[eye(r) zeros(r,Z-r)]/W*model.d;
gN_D2TUC = -(1/model.C)*(model.BG'*model.BG)\model.BG'*model.d;

%% Controller gain synthesis 
% Compute LQR weight matrices as detailled in [1] for DTUC
Q_DTUC = [eye(r) zeros(r,Z-r)]*W'*diag(1./model.capacity)*...
    W*[eye(r) ;zeros(Z-r,r)];
R_DTUC = 0.0001*eye(model.stg);

% Compute LQR weight matrices as detailled in [1] for D2TUC
Q_D2TUC = diag(1./model.capacity);
R_D2TUC = 0.0001*eye(model.L);

% Centralized gain computation for DTUC
[K_TUC,P_TUC] = LQROneStepLTI_augmented(A1_hat,Bg1_hat,Q_DTUC,R_DTUC,ones(size(model.E_DTUC_phi)),1e3,1e-5,model.A,model.Bg,r,Z,W);    
% One-step gain computation for DTUC with configuration psi
[K_DTUC_psi,P_DTUC_psi] = LQROneStepLTI_augmented(A1_hat,Bg1_hat,Q_DTUC,R_DTUC,model.E_DTUC_psi,1e3,1e-5,model.A,model.Bg,r,Z,W);
% One-step gain computation for DTUC with configuration phi
[K_DTUC_phi,P_DTUC_phi] = LQROneStepLTI_augmented(A1_hat,Bg1_hat,Q_DTUC,R_DTUC,model.E_DTUC_phi,1e3,1e-5,model.A,model.Bg,r,Z,W);
% Centralized gain computation for D2TUC with configuration phi
[K_D2TUC_C,P_D2TUC_C] = LQRCentralizedLTI(model.A,model.BG,Q_D2TUC,R_D2TUC);
% One-step gain computation for D2TUC with configuration psi
[K_D2TUC_psi,P_D2TUC_psi] = LQROneStepLTI(model.A,model.BG,Q_D2TUC,R_D2TUC,model.E_D2TUC_psi);
% One-step gain computation for D2TUC with configuration phi
[K_D2TUC_phi,P_D2TUC_phi] = LQROneStepLTI(model.A,model.BG,Q_D2TUC,R_D2TUC,model.E_D2TUC_phi);
 
%% Nonlinear simulation
controlStrat = 6; % Number of control strategies to simulate
nDisc = 10; % Number of discrete time steps to simulate
tspan = (0:model.T:nDisc*model.C); % (s) 
xNL = cell(controlStrat,1);
xDisc = cell(controlStrat,1);
dNL = cell(controlStrat,1); 
gNL = cell(controlStrat,1);
uNL = cell(controlStrat,1);

% Variable initialization
for m = 1:controlStrat 
    xNL{m,1} = zeros(model.L,length(tspan));
    dNL{m,1} = zeros(model.L,length(tspan)-1);
    gNL{m,1} = zeros(model.stg,length(0:model.C:tspan(end))-1);
    uNL{m,1} = zeros(model.stg,length(tspan(end))-1);
end

% Simulate each control strategy
for  m = 1:controlStrat
    xNL{m,1}(:,1) = model.x0;
    for k = 1:length(tspan)-1
        % Control update frequency is T/C times slower if
        if ~rem(int16(k-1),int16(model.C/model.T))
            if k ~= 1
               xDisc{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) = xNL{m,1}(:,k);
            else
                xDisc{m,1}(:,1) = model.x0;
            end
            xD = xDisc{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1);
            switch m 
            case 1 % TUC
                gNL{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) =...
                    LQcontrolAction(xD,K_TUC,model,gN_DTUC);
            case 2 % DTUC with configuration psi     
                gNL{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) =...
                    LQcontrolAction(xD,K_DTUC_psi,model,gN_DTUC);
            case 3 % DTUC with configuration phi                 
                gNL{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) =...
                    LQcontrolAction(xD,K_DTUC_phi,model,gN_DTUC);
            case 4 % TUC - QPC (centralized version of D2TUC) inspired in [2] 
                gNL{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) =...
                    QPCcontrolAction(xD,K_D2TUC_C,model,model.stageMatrix,gN_D2TUC);
            case 5 % DTUC with configuration phi  
                gNL{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) =...
                    QPCcontrolAction(xD,K_D2TUC_psi,model,model.stageMatrix,gN_D2TUC);
            case 6 % DTUC with configuration phi
                gNL{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) =...
                    QPCcontrolAction(xD,K_D2TUC_phi,model,model.stageMatrix,gN_D2TUC);
            end
        end

        % Compute nonlinear control action with upstream gating
        for l = 1:model.L % Equation (14) of [2]
            if sum(xNL{m,1}(model.turningRatesTable(l,1:end-1)~=0,k) >=...
                   model.jam*model.capacity(model.turningRatesTable(l,1:end-1)~=0))
               uNL{m,1}(l,k) = 0;
            else
               uNL{m,1}(l,k) = min(xNL{m,1}(l,k)/model.T,model.stageMatrix(l,:)*...
                   gNL{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1)*...
                model.saturation(l)/model.C);
            end
        end
        dNL{m,1}(:,k) = min((model.capacity-...
            xNL{m,1}(:,k)-model.Bu_sim*uNL{m,1}(:,k))/model.T,model.d);
        xNL{m,1}(:,k+1) = xNL{m,1}(:,k)+model.Bu_sim*uNL{m,1}(:,k)+model.T*model.d;
        
        % Catch overspill
        if sum(xNL{m,1}(:,k+1)./model.capacity>1) ~=0
           fprintf("Overspill: method %d | instant k=%d | link %d \n", m,k+1,find(xNL{m,1}(:,k+1)./model.capacity>1));
           break;
        end
      
    end
end

%% Compute performance indices defined in [1,2]
TTS = zeros(controlStrat,1);
RBQ = zeros(controlStrat,1);
xNL_crit = cell(controlStrat,1);
for  m = 1:controlStrat
    xNL_crit{m,1}(:,1) = model.x0;
    for k = 2:length(tspan)-1
        % Control update frequency is T/C times slower if
        if ~rem(int16(k-1),int16(model.C/model.T))
            xNL_crit{m,1}(:,idivide(int16(k-1),int16(model.C/model.T))+1) = ...
                 mean(xNL{m,1}(:,k-(model.C/model.T-1):k),2);
        end
    end
end
for m = 1:controlStrat % variable initialization
    TTS(m) = model.C*(1/3600)*sum(sum(xNL_crit{m,1}(:,2:end)));
    RBQ(m) = sum(sum(xNL_crit{m,1}(:,2:end).^2,2)./model.capacity);
end

%% Comparison
tb = [TTS';zeros(1,6);RBQ';zeros(1,6)];
for i = [2 3 5 6]
   if (i==3) || (i==2)
        tb(2,i) = 100*(tb(1,i)-tb(1,1))/tb(1,1);
        tb(4,i) = 100*(tb(3,i)-tb(3,1))/tb(3,1);
   else
       tb(2,i) = 100*(tb(1,i)-tb(1,4))/tb(1,4);
       tb(4,i) = 100*(tb(3,i)-tb(3,4))/tb(3,4);
   end
end
fprintf("--------------------------------------------------------------------------------------------\n")
for i = 1:4 
    if rem(i,2)==1
        fprintf("%g\t%.4g\t%.4g\t%g\t%.4g\t%.4g\n",tb(i,1),tb(i,2),tb(i,3),tb(i,4),tb(i,5),tb(i,6));
    else
        fprintf("--\t%.3g\t%.3g\t--\t%.3g\t%.3g\n",tb(i,2),tb(i,3),tb(i,5),tb(i,6));
        fprintf("--------------------------------------------------------------------------------------------\n")
    end 
end

%% Plots
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
for m = 1:6
    p = plot(tspan/model.C,sum(abs((V*V')*xNL{m,1}(:,:))));
    p.LineWidth = 3;
end
legend('TUC', 'DTUC|\Psi', 'DTUC|\Phi','D2TUC|Cent.','D2TUC|\Psi','D2TUC|\Phi');
%set(gcf, 'Position', [100 100 900 550]);
ylabel("$\sum_i|[\mathbf{z_2}]_i|$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
hold off;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
for m = 1:6
    p = plot(tspan/model.C,sum(abs((H*H')*xNL{m,1}(:,:))));
    p.LineWidth = 3;
end
legend('TUC', 'DTUC|\Psi', 'DTUC|\Phi','D2TUC|Cent.','D2TUC|\Psi','D2TUC|\Phi');
%set(gcf, 'Position', [100 100 900 550]);
ylabel("$\sum_i|[\mathbf{z_1}]_i|$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
hold off;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
for i = 1:11
    r = plot(tspan/model.C,xNL{3,1}(i,:)./model.capacity(i),'LineWidth',3);
end
legend('z=1', 'z=2', 'z=3','z=4','z=5','z=6','z=7','z=8','z=9','z=10','z=11');
%set(gcf, 'Position', [100 100 900 550]);
ylabel("$x_z/x_{z,max}$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
hold off;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
for i = 1:9 %[5 8 9]
    r = plot((0:model.T/100:nDisc*model.C-model.T/2)/model.C,gNL{3,1}(i,idivide(int16(((model.T/2:model.T/100:nDisc*model.C))/model.T-1),int16(model.C/model.T))+1),'LineWidth',3);
    %r.LineWidth = 3;
end
legend('s=1', 's=2', 's=3','s=4', 's=5', 's=6','s=7', 's=8', 's=9');
%set(gcf, 'Position', [100 100 900 550]);
ylabel("$g_s$ (s)",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
hold off;



%% Auxiliary functions

%% modelSynthesis - Description
% This function synthesizes the model of a traffic network from the raw
% data provide in text files. 
% Input:    - folder: path to directory of text files
%               - 'general.txt': general parameters
%               - 'junctions_table.txt': caracteristics of each junction
%               - 'links_table.txt': characteristics of each link
%               - 'stage_matrix.txt': stage matrix
%               - 'turning_rates_table.txt': turning rates matrix
% (For templates of the text files check "template_*.txt")
% Output:   - model: struct of variables that characterize the network
function model =  modelSynthesis(folder)
%% Import general parameters
data = importdata(folder+"general.txt");
J = data(1,1); % number of junctions (nodes)
L = data(1,2); % number of links
stg = data(1,3); % number of stages
C = data(1,4); % control cycle time-step (s)
H = data(1,5);
jam = data(1,6); % factor c. threshold for definition of a jam
T = data(1,7); % simulation time-step (s)

%% Import junctions table
junctionTable = importdata(folder+"junctions_table.txt");
lostTime = junctionTable(:,1); % time lost at each junction during a control cycle (s)
stages = junctionTable(:,2); % number of stages of each junction
gmin = junctionTable(:,3); % minimum green time for each stage at each junction (s)

%% Import links table
linksTable = importdata(folder+"links_table.txt");
capacity = linksTable(:,1); % maximum nuber of vehicles in each link
saturation = linksTable(:,2)/(60^2); % saturation flows (ve/s)
lanes = linksTable(:,3); % number of lanes of each link
x0 = linksTable(:,4); % initial condition
d = linksTable(:,5)/(60^2); % demand flows (ve/s)

%% Import stage matrix
stageMatrix = importdata(folder+"stage_matrix.txt"); % stage matrix (ROW) 

%% Import turning rates matrix
turningRatesTable = importdata(folder+"turning_rates_table.txt"); 

%% Traffic network graph variables
% Cell to hold the stages of each junction
junctions = cell(J,1); 
count = 1;
% Matrix to hold the junction of origin and destination of each link
links = zeros(L,2);
% Fill variables 'junctions' and 'links'
for j = 1:J
    junctions{j,1} = count:count+stages(j)-1; % stages of junction
    count = count+stages(j);
end
for l = 1:L
    stageOfLink = find(stageMatrix(l,:)==1);
    for j = 1:J
        if sum(intersect(junctions{j,1},stageOfLink))
            links(l,2) = j;
            break;
        end
    end
end
for l = 1:L
    stageOfLink = find(stageMatrix(l,:)==1);
    for j = 1:J
        if sum(intersect(junctions{j,1},stageOfLink))
            links(l,2) = j;
            break;
        end
    end
end
for l = 1:L
    linkFrom = find(turningRatesTable(:,l)~=0);
    for lF = 1:length(linkFrom)
        links(l,1) = links(linkFrom(lF),2);
    end
end
% Variable to hold the indices of the links that originate outside the network
inLinks = find(links(:,1)==0);
% Variable to hold the indices of the links that do not originate outside the network
notInLinks = find(links(:,1)~=0);

%% Compute LTI state-space system matrices
% Computation according to [1]
A = eye(L);
% Linear simulation model / control model
% T_sim
Bu_sim = T*(diag(1-turningRatesTable(:,end))*turningRatesTable(:,1:end-1)'-eye(L)); 
BG_sim = (1/C)*Bu_sim*diag(saturation); 
Bg_sim = BG_sim*stageMatrix; 
% T_ctrl = C
Bu = C*(diag(1-turningRatesTable(:,end))*turningRatesTable(:,1:end-1)'-eye(L)); 
BG = (1/C)*Bu*diag(saturation); 
Bg = BG*stageMatrix; 
C_z = eye(L);

%% Compute sparcity patterns
% Sparcity pattern for method DTUC and configuration Psi
E_DTUC_psi = zeros(stg,L);
% Sparcity pattern for method DTUC and configuration Phi
E_DTUC_phi = zeros(stg,L);

% Compute sparcity patterns 
for stage = 1:stg % for all stages
    junctionOfStage = 0; % find the corresponding junction
    for j = 1:J
       if sum(junctions{j,1}==stage)
           junctionOfStage = j; % juncrion of stage s found
           break;
       end
    end
    linksFromJunction = find(links(:,1)==junctionOfStage); % find the links from juntion of stage
    E_DTUC_psi(stage,links(:,1)==junctionOfStage) = 1;
    linksTowardsJunction = find(links(:,2)==junctionOfStage); % find the links from juntion of stage
    E_DTUC_psi(stage,links(:,2)==junctionOfStage) = 1;
    E_DTUC_phi(stage,:) = E_DTUC_psi(stage,:);
    junctionOfLinkTowards = links(linksTowardsJunction,1);
    junctionOfLinkTowards = junctionOfLinkTowards(junctionOfLinkTowards~=0);
    for junctionTowardsIdx = 1:length(junctionOfLinkTowards)
        E_DTUC_phi(stage,links(:,2)==junctionOfLinkTowards(junctionTowardsIdx)) = 1;
        E_DTUC_phi(stage,links(:,1)==junctionOfLinkTowards(junctionTowardsIdx)) = 1;
    end
    junctionOfLinkFrom = links(linksFromJunction,2);
    for junctionFromIdx = 1:length(junctionOfLinkFrom)
        E_DTUC_phi(stage,links(:,2)==junctionOfLinkFrom(junctionFromIdx)) = 1;
        E_DTUC_phi(stage,links(:,1)==junctionOfLinkFrom(junctionFromIdx)) = 1;
    end

end

% Sparcity pattern for method D2TUC and configuration Psi
E_D2TUC_psi = zeros(L,L);
% Sparcity pattern for method D2TUC and configuration Phi
E_D2TUC_phi = zeros(L,L);

% Compute sparcity patterns
for link = 1:L
    junctionOfLink = links(link,2);  
    linksFromJunction = find(links(:,1)==junctionOfLink); % find the links from juntion of stage
    E_D2TUC_psi(link,links(:,1)==junctionOfLink) = 1;
    linksTowardsJunction = find(links(:,2)==junctionOfLink); % find the links from juntion of stage
    E_D2TUC_psi(link,linksTowardsJunction) = 1;
    E_D2TUC_phi(link,:) = E_D2TUC_psi(link,:);
    junctionOfLinkTowards = links(linksTowardsJunction,1);
    junctionOfLinkTowards = junctionOfLinkTowards(junctionOfLinkTowards~=0);
    for junctionTowardsIdx = 1:length(junctionOfLinkTowards)
        E_D2TUC_phi(link,links(:,2)==junctionOfLinkTowards(junctionTowardsIdx)) = 1;
        E_D2TUC_phi(link,links(:,1)==junctionOfLinkTowards(junctionTowardsIdx)) = 1;
    end
    junctionOfLinkFrom = links(linksFromJunction,2);
    for junctionFromIdx = 1:length(junctionOfLinkFrom)
        E_D2TUC_phi(link,links(:,2)==junctionOfLinkFrom(junctionFromIdx)) = 1;
        E_D2TUC_phi(link,links(:,1)==junctionOfLinkFrom(junctionFromIdx)) = 1;
    end
end

%% Fill model struct
model = struct('junctionTable',junctionTable,...
               'linksTable',linksTable,...
               'stageMatrix',stageMatrix,...
               'turningRatesTable',turningRatesTable,...
               'J',J,...
               'L',L,...
               'stg',stg,...
               'C',C,...
               'H',H,...
               'jam',jam,...
               'T',T,...
               'lostTime',lostTime,...
               'stages',stages,...
               'gmin',gmin,...
               'capacity',capacity,...
               'saturation',saturation,...
               'lanes',lanes,...
               'x0',x0,...
               'd',d,...%'junctions',junctions,...               
               'links',links,...
               'inLinks',inLinks,...
               'notInLinks',notInLinks,...
               'A',A,...
               'Bu',Bu,...
               'BG',BG,...
               'Bg',Bg,...
               'Bu_sim',Bu_sim,...
               'BG_sim',BG_sim,...
               'Bg_sim',Bg_sim,...
               'C_z',C_z,...
               'E_DTUC_psi',E_DTUC_psi,...
               'E_DTUC_phi',E_DTUC_phi,...
               'E_D2TUC_psi',E_D2TUC_psi,...
               'E_D2TUC_phi',E_D2TUC_phi);
end

%% LQcontrolAction - Description
% This function implements a standard linear quadratic feedback control
% action for the raffic network model
% Input:    - x: state
%           - L_LQ: LQ gain matrix
%           - model: struct of variables that characterize the network
%           - gN: historic green times
% Output:   - g: green times for each stage
function g = LQcontrolAction(x,L_LQ,model,gN)
    gpre = -L_LQ*x+gN;
    g = zeros(size(gpre));
    for j = 1:model.J   
        stgsj = (sum(model.stages(1:j-1))+1:sum(model.stages(1:j-1))+model.stages(j));
        a = (gpre(stgsj)-model.gmin(j));
        d = ones(model.stages(j),1);
        c = model.C-model.lostTime(j)-model.stages(j)*model.gmin(j);
        b = c*ones(model.stages(j),1);
        g(stgsj) = knapsack(a,b,c,d)+model.gmin(j);     
    end
end

%% QPCcontrolAction - Description
% This function implements a linear quadratic rogramming feedback control
% action for the raffic network model
% Inspired in the nonlinear QPC control action proposed in [2]
% Input:    - x: state
%           - L_QPC: QPC linear gain matrix
%           - model: struct of variables that characterize the network
%           - ROW: stage matrix S
%           - gNQPC: historic green times
% Output:   - g: green times for each stage
function g = QPCcontrolAction(x,L_QPC,model,ROW,gNQPC)
    Gpre = -L_QPC*x+gNQPC;
    gpre = (ROW'*ROW)\ROW'*Gpre;
    g = zeros(size(gpre));
    for j = 1:model.J   
        stgsj = (sum(model.stages(1:j-1))+1:sum(model.stages(1:j-1))+model.stages(j));
        a = (gpre(stgsj)-model.gmin(j));
        d = ones(model.stages(j),1);
        c = model.C-model.lostTime(j)-model.stages(j)*model.gmin(j);
        b = c*ones(model.stages(j),1);
        g(stgsj) = knapsack(a,b,c,d)+model.gmin(j);     
    end
end

%% knapsack - Description
% This function solves the knapsack problem 
% Implementation of algorithm in [3]
% Input:    - a,b,c,d as defined in [1]
% Output:   - x: knapsack solution
function x = knapsack(a,b,c,d)
% Variables
aux = a-d.*b;
aux = [aux;a];
y = sort(aux);
n = length(a);
% 0. Initialization
if c<0 || c>sum(b)
    x = NaN;
    return;
else
    l = 1;
    r = 2*n;
    R = 0;
    L = sum(b);
end
% 1. Test for bracketing
while true
    if r-l == 1
        % 5. Interpolate
        lambda = y(l) + (y(r)-y(l))*(c-L)/(R-L);
        break;
    else
        m = floor((l+r)/2);
    end
    % 3. Compute new value
    aux = (a-y(m))./d;
    for i = 1:n
       aux(i) = max(min(aux(i),b(i)),0);
    end
    C = sum(aux);
    % 4. Update
    if C == c
        lambda = y(m);
        break;
    elseif C > c
        l = m;
        L = C;
    else
        r = m;
        R = C;
    end    
end
x = zeros(n,1);
for i = 1:n
    x(i) = max(min((a(i)-lambda)/d(i),b(i)),0);
end
end

%% LQROneStepLTI_augmented - Description
% This function computes the steady-state augmented one-step LQR regulator 
% gain for a window w. Method derived in [1].
% Input:    - A_hat, B_hat
%           - Q, R
%           - E: sparsity pattern
%           - itMax: maximum number of iterations until convergence 
%           - epslInf: minimum relative improvement
%           - A, B
%           - r,Z,W: as defined in [1]
% Output:   - K: nxo steady-state gain matrix
%           - P: nxn steady-state estimation error covariance matrix
% Important notes: 
%           - output gain corresponds to the control law: u(k)=-K(k)*x(k)
% WARNING: Returns Kinf = NaN and Pinf = NaN if convergence could not be reached 
function [K,P] = LQROneStepLTI_augmented(A_hat,B_hat,Q,R,E,itMax,epslInf,A,B,r,Z,W)
% Gain computation
n = size(E,2); % Get value of n from the size of A 
m = size(E,1); % Get value of n from the size of B  
P = Q; % terminal condition
Pprev = NaN;
it = itMax;
while it > 0 % LQ iterations
    K = zeros(m,n);
    S = R+B_hat'*(P)*B_hat;
    for i = 1:n
        L = zeros(n);
        L (i,i) = 1; % Generate matrix L_i
        M = zeros(m);
        for j = 1:m % Gererate matrix M_i
            if E(j,i) ~= 0
                M(j,j) = 1;
            end
        end
        % Compute the ith term of the summation 
        K = K + (eye(m)-M+M*S*M)\(M*(B_hat')*P*A_hat*([eye(r) zeros(r,Z-r)]/W)*L');
    end
    % Update P
    P_ =((W')\[eye(r);zeros(Z-r,r)]) *(P)* ([eye(r) zeros(r,Z-r)]/W);
    Q_ = ((W')\[eye(r);zeros(Z-r,r)]) *(Q)* ([eye(r) zeros(r,Z-r)]/W);
    P = Q_+K'*R*K+...
        (A-B*K)'*P_*(A-B*K);
    P = [eye(r) zeros(r,Z-r)]*W'*P*W*[eye(r);zeros(Z-r,r)];
    % Check convergence
    it = it-1;
    if abs(trace(P)-trace(Pprev))/trace(Pprev) < epslInf
        break; 
    end
    Pprev = P;
    if it == 0
        fprintf("One-step did not converge.\n");
        P = NaN;
        K = NaN;
    end
end 
end

%% References
% [1] Pedroso, L. and Batista, P., 2021. Decentralized store-and-forward 
% based strategies for the signal control problem in large-scale congested 
% urban road networks. Transportation Research Part C: Emerging 
% Technologies, 132, p.103412. doi:10.1016/j.trc.2021.103412.

% [2] Aboudolas, K., Papageorgiou, M., Kosmatopoulos, E., 2009. 
% Store-and-forward based methods for the signal control problem in 
% large-scale congested urban road networks. Transp. Res. C 17 (2), 163?174.

% [3] Helgason, R., Kennington, J., Lall, H., 1980. A polynomially bounded 
% algorithm for a singly constrained quadratic program. Math. Program. 
% 18 (1), 338?343.
