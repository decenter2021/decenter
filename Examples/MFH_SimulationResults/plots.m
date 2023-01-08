
%% Plot performnace small network
Wmax = 10;
figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
load('./models/model1/synth_C.mat','PC');
plot([1 Wmax],[trace(PC) trace(PC)],'--','LineWidth',3);
load('./models/model1/synth_OS.mat','POS');
plot([1 Wmax],[trace(POS) trace(POS)],'--','LineWidth',3);
load('./models/model1/synth_FH.mat','PFH');
plot([1 Wmax],[trace(PFH) trace(PFH)],'--','LineWidth',3);
tr_mfh = zeros(Wmax,1);
for W = 1:Wmax
    load(sprintf('./models/model1/synth_MFH_%d.mat',W),'PMFH');
    tr_mfh(W) = trace(PMFH);
end
plot(1:Wmax,tr_mfh,'*-','LineWidth',3,'MarkerSize',10);
ylabel("$\mathrm{tr}(\mathbf{P}_{\infty})$",'Interpreter','latex');
xlabel('$W_{ss}$','Interpreter','latex');
xlim([1 Wmax]);
ylim([0 60]);
legend({'C','OS','FH', 'MFH'},'Interpreter','latex','Location','southeast');
hold off;

%% Plot ev trace - small network

Wss = 5;
load('./models/model1/synth_MFH_5.mat','PMFH_seq')

figure;
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
trPseq = zeros(Wss+1,1);
for i = 1:Wss+1
    trPseq(i) = trace(PMFH_seq{i,1});
end
p = plot(0:Wss,trPseq,'-*');
p.LineWidth = 3;
p.MarkerSize = 15;
set(gca, 'YScale', 'log');
ylabel("$\mathrm{tr}(\mathbf{P}(\tau|\tau|"+sprintf("k=%d))$", 33),'Interpreter','latex');
xlabel('$\tau-(k-W_{ss})$','Interpreter','latex');
hold off;

%% Plot MC small network 

figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
load('./models/model1/mc.mat','tr_P_mc','methods');
plot(0:size(tr_P_mc,2)-1,tr_P_mc(methods,:)','LineWidth',3);
ylabel("$\mathrm{tr}(\mathbf{P}_{\mathrm{MC}})$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend({'C','OS','FH','MFH'},'Interpreter','latex','Location','southeast');
hold off;

perf = mean(tr_P_mc(methods,end-24:end)');
perf_rel = (perf-perf(4))/perf(4);

%% Plot performance large network
Wmax = 7;
figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
load('./models/model40/synth_C.mat','PC');
plot([1 Wmax],[trace(PC) trace(PC)],'--','LineWidth',3);
load('./models/model40/synth_OS.mat','POS');
plot([1 Wmax],[trace(POS) trace(POS)],'--','LineWidth',3);
load('./models/model40/synth_FH.mat','PFH');
plot([1 Wmax],[trace(PFH) trace(PFH)],'--','LineWidth',3);
tr_mfh = zeros(Wmax,1);
for W = 1:Wmax
    load(sprintf('./models/model40/synth_MFH_%d.mat',W),'PMFH');
    tr_mfh(W) = trace(PMFH);
end
plot(1:Wmax,tr_mfh,'*-','LineWidth',3,'MarkerSize',10);
ylabel("$\mathrm{tr}(\mathbf{P}_{\infty})$",'Interpreter','latex');
xlabel('$W_{ss}$','Interpreter','latex');
xlim([1 Wmax]);

legend({'C','OS','FH', 'MFH'},'Interpreter','latex','Location','southwest');
hold off;

%% Plot ev trace - small network

Wss = 5;
load('./models/model40/synth_MFH_5.mat','PMFH_seq')

figure;
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
trPseq = zeros(Wss+1,1);
for i = 1:Wss+1
    trPseq(i) = trace(PMFH_seq{i,1});
end
p = plot(0:Wss,trPseq,'-*');
p.LineWidth = 3;
p.MarkerSize = 15;
set(gca, 'YScale', 'log');
ylabel("$\mathrm{tr}(\mathbf{P}(\tau|\tau|"+sprintf("k=%d))$", 59),'Interpreter','latex');
xlabel('$\tau-(k-W_{ss})$','Interpreter','latex');
hold off;


%% Plot simulation uncorrelated version of model40 to see divergence of PMHE1

sim = load('./models/model41/simulation.mat');
error_norm =  sim.error_norm(sim.methods,:);

figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
p = plot(0:100,error_norm,'LineWidth',3);
ylabel("$||\hat{\mathbf{x}}(k|k|k)-\mathbf{x}(k)||_2$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend({'C','OS', 'PMHE1 ($W_{\mathrm{PMHE1}} = 5$)'},'Interpreter','latex','Location','best');
hold off;

mean_error_norm =  mean(error_norm');
mean_error_norm_rel = (mean_error_norm-mean_error_norm(3))/mean_error_norm(3);

%% Plot MC large network 

figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
load('./models/model40/mc.mat','tr_P_mc','methods');
plot(0:size(tr_P_mc,2)-1,tr_P_mc(methods,:)','LineWidth',3);
ylabel("$\mathrm{tr}(\mathbf{P}_{\mathrm{MC}})$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend({'C','OS','FH','MFH'},'Interpreter','latex','Location','southeast');
hold off;

perf = mean(tr_P_mc(methods,end-19:end)');
perf_rel = (perf-perf(4))/perf(4);



%% Simulation large weak couplings network

simN2 = load('./models/model42/simulation_N2.mat');
simN5 = load('./models/model42/simulation_N5.mat');

error_norm = zeros(5,301);
error_norm(1:2,:) = simN2.error_norm(1:2,:);
error_norm(3,:) = simN2.error_norm(4,:);
error_norm(4,:) = simN2.error_norm(5,:);
error_norm(5,:) = simN5.error_norm(5,:);

figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
p = plot(0:300,error_norm,'LineWidth',3);
ylabel("$||\hat{\mathbf{x}}(k|k|k)-\mathbf{x}(k)||_2$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend({'C','OS','MFH ($W_{ss} = 2$)', 'PMHE1 ($W_{\mathrm{PMHE1}} = 2$)', 'PMHE1 ($W_{\mathrm{PMHE1}} = 5$)'},'Interpreter','latex');
hold off;

mean_error_norm =  mean(error_norm');
mean_error_normel = (mean_error_norm-mean_error_norm(3))/mean_error_norm(3);



%% Show consistent


sim = load('./models/model40/simulation.mat');
load('models/model40/projected_ss.mat');

component = round(rand()*1000);
figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
for m = 4% sim.methods
    p = plot(0:1000,sim.error{m,1}(component,:),'LineWidth',3);
end
plot([0 1000],3*sqrt(P_ss{4,1}(component,component))*[1 1],'--','Color','k','LineWidth',3);
plot([0 1000],-3*sqrt(P_ss{4,1}(component,component))*[1 1],'--','Color','k','LineWidth',3);
ylabel(sprintf("$[\\hat{\\mathbf{x}}(k|k|k)]_{%d}-[\\mathbf{x}(k)]_{%d}$",component,component),'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend({'MFH ($W_{ss} = 5$)','$3\sigma$ bound'},'Interpreter','latex');
hold off;





%% Plots 



%% 

load('./models/model40/mc.mat');
load('./models/model40/synth_C.mat','PC');
figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
p = plot(0:200,tr_P_mc(methods,:)','LineWidth',3);
ylabel("$\mathrm{tr}(\mathbf{P}_{\mathrm{MC}}(k))$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend({'C','OS','FH','MFH ($W_{ss} = 5$)'},'Interpreter','latex','Location','best');
hold off;

perf_norm_mc = mean(tr_P_mc(methods,end-19:end),2)'/trace(PC)
(perf_norm_mc-perf_norm_mc(4))/perf_norm_mc(4)

