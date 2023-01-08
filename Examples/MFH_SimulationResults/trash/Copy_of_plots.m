
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
load('./models/model1/synth_H2.mat','P');
plot([1 Wmax],[trace(P) trace(P)],'--','LineWidth',3);
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
legend({'C','OS','FH', 'H2', 'MFH'},'Interpreter','latex','Location','southeast');
hold off;

%% Plot MC small network 

figure;
hold on;
set(gca,'FontSize',30);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
load('./models/model1/mc.mat','tr_P_mc');
plot(0:length(tr_P_mc)-1,tr_P_mc,'LineWidth',3);
ylabel("$\mathrm{tr}(\mathbf{P}_{\mathrm{MC}})$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
% xlim([1 Wmax]);
% ylim([0 60]);
legend({'C','OS','FH', 'H2', 'MFH'},'Interpreter','latex','Location','southeast');
hold off;
