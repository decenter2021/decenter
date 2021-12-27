%% Example of LQROneStepLTI
% Method developed in [1]
%% Synthetic random system
n = 5;
m = 3;
A = [0.6085    0.0188    0.9615    0.6161    0.0494;
     0.3959    0.6575    0.7305    0.5338    0.4611;
     0.1743    0.9322    0.3419    0.2229    0.1466;
     0.5204    0.7850    0.2257    0.4315    0.8538;
     0.8603    0.8842    0.1132    0.7550    0.3850];
B = [0.3843    0.7494    0.4509;
     0.1446    0.5369    0.0092;
     0.6133    0.6413    0.6322;
     0.3401    0.4020    0.3725;
     0.7084    0.4744    0.0031];
Q = [0.8057    0.9316    0.9227    0.7569    0.6049;
     0.9316    2.6200    1.2024    1.4863    1.6318;
     0.9227    1.2024    1.1229    0.9413    0.6800;
     0.7569    1.4863    0.9413    1.3068    1.2506;
     0.6049    1.6318    0.6800    1.2506    1.5626];
R = [1.0134    0.5867    0.9654;
     0.5867    0.4666    0.4427;
     0.9654    0.4427    1.5383];
E = [1     1     1     0     1;
     0     0     1     0     0;
     0     0     1     1     1];

%% Gain synthesis 
opts.verbose = true;
[Kinf,Pinf] = LQROneStepLTI(A,B,Q,R,E,opts);
Kinf
trace(Pinf)
 
%% Simulate error dynamics
% Generate random initial covariance
P0 = rand()*eye(n);
% Simulation time
SimIt = 50;
% Initialise error cell
x = cell(1,SimIt);
% Generate random initial error
x0 = transpose(mvnrnd(zeros(n,1),P0));
for j = 1:SimIt
    if j == 1
        x{1,j} = (A-B*Kinf)*x0;
    else
        x{1,j} = (A-B*Kinf)*x{1,j-1};
    end
end

%% Plot dynamics
% Plot the ||x||_2 vs instant of the simulation 
figure;
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
xPlot = zeros(SimIt,1);
for j = 1:SimIt
   xPlot(j,1) =norm(x{1,j}(:,1)); 
end
plot(0:SimIt, [norm(x0(:)); xPlot(:,1)],'LineWidth',3);
set(gcf, 'Position', [100 100 900 550]);
ylabel('$\|\mathbf{x}_{OS}(k)\|_2$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
hold off;

%% References
% [1] Viegas D, Batista P, Oliveira P, Silvestre C. Distributed controller design 
% and performance optimization for discrete-time linear systems. Optim Control 
% Appl Meth. 2020;1-18. https://doi.org/10.1002/oca.2669