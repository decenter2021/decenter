
% Wss = 5;
% [K,Pinf,Pseq] = synthesis_mfh(0,1e-4,Wss,100);
% save("./ev_trace_mfh.mat",'K','P','Pinf');


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