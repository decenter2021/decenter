%% Plots for DEKFConstellation Cart
%% Init
sat = 1;
%% Plot cartesian error with 3 sigma bounds
T = size(error{sat,1},2)-1;
figure('Position',[0 0 1000 500]);
hold on;
grid on; 
set(gca,'FontSize',15)
plot(0:T-1,error{sat,1}(1,1:end-1),'LineWidth',3);
plot(0:T-1,error{sat,1}(2,1:end-1),'LineWidth',3);
plot(0:T-1,error{sat,1}(3,1:end-1),'LineWidth',3);
plot(0:T-1,3*sqrt(max(P_pos_log(sat,1:end,:),[],3)),'--','LineWidth',3,'Color','k')
plot(0:T-1,-3*sqrt(max(P_pos_log(sat,1:end,:),[],3)),'--','LineWidth',3,'Color','k')
xlabel('$t (\mathrm{s})$','Interpreter','latex');
ylabel(sprintf('$\\hat{\\mathbf{p}}(t)-\\mathbf{p}(t) \\;(\\mathrm{m})\\;\\;$ [Sattelite \\#%05d]',sat),'Interpreter','latex')
legend({'$\mathbf{o_x}$','$\mathbf{o_y}$','$\mathbf{o_z}$','Maximum $3\sigma$-bound'},'Interpreter','latex')
hold off

%% Plot trace vs number of satellites within ISL range 

no_sats = ones(T,1);
for i = 2:T
    no_sats(i) = size(FimHist{i,1}{sat,1},1)-1;
end

figure('Position',[0 0 1000 500]);
hold on;
grid on; 
set(gca,'FontSize',15)
xlabel('$t (\mathrm{s})$','Interpreter','latex');

yyaxis left;
plot(0:T-1,trace_log(sat,1:end),'LineWidth',3)
ylabel(sprintf('$\\mathrm{tr}\\left(\\mathbf{P}_{i,(i,i)}(t)\\right)\\;$ [Sattelite \\#%05d]',sat),'Interpreter','latex')
yyaxis right
plot(0:T-1,no_sats,'LineWidth',3)
ylabel(sprintf('Number of satellites within ISL range'),'Interpreter','latex')
hold off

%% RMSE 
sat = 1;
T = size(error{sat,1},2)-1;
W = 301;
rmse = zeros(T,1);
for t = 1:T
    for i = max([1 t-(W-1)/2]):min([T t+(W-1)/2])
        rmse(t) = rmse(t) + error{sat,1}(1:3,i)'*error{sat,1}(1:3,i);    
    end
    rmse(t) = sqrt((1/(min([T t+(W-1)/2])-max([1 t-(W-1)/2])+1))*rmse(t));
end

figure('Position',[0 0 1000 500]);
hold on;
grid on; 
set(gca,'FontSize',15)
plot(0:T-1,rmse,'LineWidth',3);
xlabel('$t (\mathrm{s})$','Interpreter','latex');
ylabel(sprintf('$\\mathrm{RMSE}\\left(||\\hat{\\mathbf{p}}(t)-\\mathbf{p}(t)|| \\right)\\;(\\mathrm{m})\\;\\;$ [Sattelite \\#%05d]',sat),'Interpreter','latex')
hold off
