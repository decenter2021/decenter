%% Example of kalmanFiniteHorizonLTI
%% Synthetic system of Viegas et al. (2018) 
n = 5;
o = 4;
A = [0.152  0.092   0.235   0.642   0.506;
     0.397  0.615   0.448   0.221   0.279;
     0.375  0.011   0.569   0.837   0.747;
     0.131  0.573   0.061   0.971   0.237;
     0.435  0.790   0.496   0.846   0.957];
C = [0.620  0.255   0.725   0.404   0.511;
     0.600  0.859   0.230   1.988   0.061;
     0.173  0.911   0.576   0.090   0.726;
     0.090  0.700   0.811   0.321   0.557];
Q = [3.318  4.662   1.598   -1.542  -1.999;
     4.662  11.520  2.608   -2.093  -5.442;
     1.598  2.608   4.691   0.647   -0.410;
     -1.542 -2.093  0.647   2.968   0.803;
     -1.999 -5.442  -0.410  0.803   2.851];
R = [3.624  2.601   -0.042  -0.944;
     2.601  7.343   -0.729  -2.786;
     -0.042 -0.729  0.745   -0.242;
     -0.944 -2.786  -0.242  1.612];
E = [1   0   1   1;
     0   1   0   1;
     0   0   1   0;
     1   1   1   0;
     1   1   0   1];
 %% Gain synthesys 
opts.verbose = true;
opts.W = 30;
opts.maxOLIt = 10;
[Kinf,Pinf] = kalmanFiniteHorizonLTI(A,C,Q,R,E,opts);
Kinf 
trace(Pinf)
%% Simulate error dynamics
% Generate random initial covariance
P0 = rand(n,n);
P0 = 100*(P0*P0');
% Simulation time
SimIt = 100;
% Initialise error cell
error = cell(1,SimIt);
% Generate random initial error
error0 = transpose(mvnrnd(zeros(n,1),P0));
for j = 1:SimIt
    if j == 1
        error{1,j} = (eye(n)-K*C)*(A*error0+...
            mvnrnd(zeros(n,1),Q)')-K*mvnrnd(zeros(o,1),R)';
    else
        error{1,j} = (eye(n)-K*C)*(A*error{1,j-1}+...
            mvnrnd(zeros(n,1),Q))'-K*mvnrnd(zeros(o,1),R)';
    end
end

%% Plot estimation error norm
% Plot the ||error||_2 vs instant of the simulation 
figure;
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
%title('Norm of estimation error simulation - centralized gain'); 
errorPlot = zeros(SimIt,1);
for j = 1:SimIt
   errorPlot(j,1) =norm(error{1,j}(:,1)); 
end
plot(0:SimIt, [norm(error0(:)); errorPlot(:,1)],'LineWidth',3);
set(gcf, 'Position', [100 100 900 550]);
ylabel('$\|\hat{\mathbf{x}}_{FH}(k|k)-\mathbf{x}(k)\|_2$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
hold off;