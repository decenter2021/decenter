%% Tutorial of LQRCentralizedLTV
%% Synthetic random system
T = 50;
n = 5;
m = 3;
rng(1); % Pseudo-random seed for consistency
% Alternatively comment out rng() to generate a random system
% Do not forget to readjust the synthesys parameters of the methods
system = cell(T+1,4);
for i = 1:T+1
    if i == 1
        system{i,1} = rand(n,n)-0.5;
        system{i,2} = rand(n,m)-0.5;
    elseif i == T+1
        system{i,1} = nan;
        system{i,2} = nan;
        system{i,3} = rand(n,n)-0.5;
        system{i,3} = system{i,3}*system{i,3}';
        system{i,4} = nan;
        continue;
    else % Generate time-varying dynamics preventing erratic behaviour
        system{i,1} = system{i-1,1}+(1/4)*(rand(n,n)-0.5);
        system{i,2} = system{i-1,2}+(1/4)*(rand(n,m)-0.5);
    end
    system{i,3} = rand(n,n)-0.5;
    system{i,3} = system{i,3}*system{i,3}';
    system{i,4} = rand(m,m)-0.5;
    system{i,4} = system{i,4}*system{i,4}';
end

%% Synthesize regulator gain using the centralized method (with some optional parameters)
opts.verbose = true;
[K,P] = LQRCentralizedLTV(system,T,opts);
 
%% Simulate error dynamics
% Generate random initial covariance for the initial time instant 
% (LQROneStepLTV does not make assumptions on P0)
P0 = rand(n,n);
P0 = P0*P0';
% Initialise error cell
x = cell(T+1,1);
% Generate random initial error
x{1,1} = transpose(mvnrnd(zeros(n,1),P0));
for j = 2:T+1
    x{j,1} = (system{j-1,1}-system{j-1,2}*K{j-1,1})*x{j-1,1};
end

%% Plot the norm of the estimation error
% Plot the ||x||_2 vs instant of the simulation 
figure;
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
xPlot = zeros(T+1,1);
for j = 1:T+1
   xPlot(j,1) =norm(x{j,1}(:,1)); 
end
plot(0:T, xPlot(:,1),'LineWidth',3);
set(gcf, 'Position', [100 100 900 550]);
ylabel('$\|\mathbf{x}_{C}(k)\|_2$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
hold off;