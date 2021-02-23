%% Description
% This script simulates the original non-linear quadruple tank. The control
% vector is computed using predefined law (since we only want to study the 
% performance estimation). In this script observers of the LTV system are 
% also implemented. It includes a centralized observer and two distributed 
% observers, whose gains are computed using the One Step and Delayed Finite Horizon Algorithms. 

%% Descentralization scheme
% The distributed scheme presented assumes the control vector is determined
% in a centralized manner and  is available to the local observers via 
% communication. It consists of 4 local observers on in each tank,
% estimating their own position

%% Constants 
rng(1);
clear;
Cte = getConstantsWL();

%% Simulation options 
% Select initial state vector in cm
x0 = 20*ones(4,1);
% Select the time window for the simulation in seconds
tspan = [0 600];
% Compute number of discrete iterations
SimIt = tspan(2)/Cte.dT+1;

%% Reference to lower tanks
% Initialize the reference vector
T_ref = 200; % times dT
ref = cell(1,SimIt);
refEq = cell(2,SimIt);
% Step references
% Gererate random references with the defined periodicity 
aux = [20; 20];
for i = 1:T_ref:SimIt
    if aux(1) == 20 aux = [30;30]; else aux = [20;20]; end
    systemRef = getDiscreteDynamicsWL([aux;0;0],Cte);
    for j = i:min(i+T_ref-1,SimIt+1)
        ref{1,j} = zeros(Cte.n,1);
        ref{1,j}(1:Cte.n/2,1) = aux(:,1); 
        refEq{1,j} = systemRef{1,5};
        refEq{2,j} = systemRef{1,6};
    end
end
j = SimIt;
ref{1,j} = ref{1,j-1};
refEq{1,j} = refEq{1,j-1};
refEq{2,j} = refEq{2,j-1};
j = SimIt+1;
ref{1,j} = ref{1,j-1};
refEq{1,j} = refEq{1,j-1};
refEq{2,j} = refEq{2,j-1};

%Sinusoidal references
%Gererate random references with the defined periodicity 
for i = 1:SimIt+1
    ref{1,i}(2,1) = 30+10*cos((i-1)/35);
    systemRef = getDiscreteDynamicsWL(ref{1,i},Cte);
    refEq{1,i} = systemRef{1,5};
    refEq{2,i} = systemRef{1,6};
end

%% Initialize controller

% Initialize temporary cell to hold the covariance and gain matrices
K = cell(4,SimIt-1);

%% Discrete simulation

% Discrete time vector
t_disc = 0:Cte.dT:tspan(2);
% Initialize discrete control vector
uControlDisc = cell(4,SimIt-1);
uEqDisc = cell(4,SimIt-1);
xEqDisc = cell(4,SimIt-1);
% Initialize measurement and state vectors
xNonLin = cell(4,SimIt);
x = cell(4,SimIt);
% Initialize nonlinear solver 
nonLinSol = cell(4,1);



for i = 1:SimIt
    i
    if i == 1     
        for m = 1:4
             % Measure output
             x{m,i}(1:Cte.n,1) = x0;
             x{m,i}(Cte.n+1:3*Cte.n/2) = x0(1:Cte.n/2,1)-ref{1,i}(1:Cte.n/2,1);
             xNonLin{m,i} = x0;
             if rem(i-1,Cte.d(m)) == 0
                [K(m,i:min(i+Cte.d(m)-1,SimIt-1)),uEqDisc(m,i:min(i+Cte.d(m)-1,SimIt-1))]...
                    = iLQR(x{m,i},refEq(:,i:end),Cte,min(Cte.d(m),SimIt-i),m);
             end
             uControlDisc{m,i} = -transpose(K{m,i})*(x{m,i}-refEq{1,i})+uEqDisc{m,i};
             uControlDisc{m,i}(uControlDisc{m,i}<0) = 0;
             uControlDisc{m,i}(uControlDisc{m,i}>Cte.uMax) = Cte.uMax;
             % Simulate nonlinear dynamics
             nonLinSol{m,1} = ode45(@(t,x) xdotContinuous(x,uControlDisc{m,min(floor(t/Cte.dT)+1,round(t_disc(i+1)/Cte.dT))},Cte),[t_disc(i) t_disc(i+1)],x0);
        end
    elseif i ~= SimIt
         for m = 1:4
             % Measure output
             x{m,i}(1:Cte.n,1) = deval(nonLinSol{m,1}, t_disc(i));
             x{m,i}(Cte.n+1:3*Cte.n/2) = x{m,i-1}(Cte.n+1:3*Cte.n/2)+x{m,i}(1:Cte.n/2,1)-ref{1,i}(1:Cte.n/2,1);
         
             % anti windup
             for j = 1:Cte.n/2
                if abs(x{m,i}(Cte.n+j)) > Cte.AntiWU(m)
                    x{m,i}(Cte.n+j) = Cte.AntiWU(m)*abs(x{m,i}(Cte.n+j))/x{m,i}(Cte.n+j); 
                end
             end
             
             xNonLin{m,i} = deval(nonLinSol{m,1}, t_disc(i));
             if rem(i-1,Cte.d(m)) == 0
                [K(m,i:min(i+Cte.d(m)-1,SimIt-1)),uEqDisc(m,i:min(i+Cte.d(m)-1,SimIt-1))]...
                    = iLQR(x{m,i},refEq(:,i:end),Cte,min(Cte.d(m),SimIt-i),m);
             end

             uControlDisc{m,i} = -transpose(K{m,i})*(x{m,i}-refEq{1,i})+uEqDisc{m,i};
             uControlDisc{m,i}(uControlDisc{m,i}<0) = 0;
             uControlDisc{m,i}(uControlDisc{m,i}>Cte.uMax) = Cte.uMax;
             % Simulate nonlinear dynamics
             nonLinSol{m,1} = ode45(@(t,x) xdotContinuous(x,uControlDisc{m,min(floor(t/Cte.dT)+1,round(t_disc(i+1)/Cte.dT))},Cte),[t_disc(i) t_disc(i+1)],deval(nonLinSol{m,1}, t_disc(i)));
         end
    else
        for m = 1:4
            xNonLin{m,i} = deval(nonLinSol{m,1}, t_disc(i));
        end
    end
end



%% Plots

dif = zeros(4,SimIt,2);

for tank = 1:4
    figure;
    xAux = zeros(4,SimIt);
    refAux = zeros(1,SimIt);
    for i = 1:SimIt
        for m = 1:4
            xAux(m,i) = xNonLin{m,i}(tank,1);
            if tank <= Cte.n/2
                refAux(1,i) = ref{1,i}(tank,1);
                dif(tank,i,m) = xAux(m,i)-ref{1,i}(tank,1);
            end
        end
    end
    hold on;
    set(gca,'FontSize',35);
    ylabel(sprintf("$h_%d$ (cm)",tank),'Interpreter','latex');
    xlabel('$t$ (s)','Interpreter','latex');
    for m = 1:4
        p = plot(t_disc,xAux(m,:));
        p.LineWidth = 3;
    end
    if tank <= Cte.n/2
         p = plot(t_disc,refAux(1,:));
         p.LineWidth = 3;
        legend(sprintf('Centralized (d = %d)',Cte.d(1)),sprintf('Centralized (d = %d)',Cte.d(2)),...
            sprintf('One-step (d = %d)',Cte.d(3)),sprintf('One-step (d = %d)',Cte.d(4)),'Reference');
    else
        legend(sprintf('Centralized (d = %d)',Cte.d(1)),sprintf('Centralized (d = %d)',Cte.d(2)),...
            sprintf('One-step (d = %d)',Cte.d(3)),sprintf('One-step (d = %d)',Cte.d(4)));
    end
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    hold off;   
    
%     if tank <= Cte.n/2
%        figure;
%        hold on;
%        for m = 1:4
%           plot(t_disc, dif(tank,:,m));
%        end
%        hold off;
%     end
end

for pump = 1:2
    figure;
    uAux = zeros(4,SimIt-1);
    for i = 1:SimIt-1
        for m = 1:4
            uAux(m,i) = uControlDisc{m,i}(pump,1);
        end
    end
    hold on;
    set(gca,'FontSize',35);
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ylabel(sprintf("$u_%d$ (V)",pump),'Interpreter','latex');
    xlabel('$t$ (s)','Interpreter','latex');
    for m = 1:4
        p = plot(t_disc(1:end-1),uAux(m,:));
        p.LineWidth = 3;
    end
    legend(sprintf('Centralized (d = %d)',Cte.d(1)),sprintf('Centralized (d = %d)',Cte.d(2)),...
            sprintf('One-step (d = %d)',Cte.d(3)),sprintf('One-step (d = %d)',Cte.d(4)));
end
