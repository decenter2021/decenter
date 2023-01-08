
%%
model = 40;
parameters = struct('sim_T',100,'MFH_w_ss',5);
methods = [true; true; true; true; false];
simulation(model, methods, parameters);

%% MC

model = 40;
parameters = struct('sim_T',100,'MFH_w_ss',5,'mc_it',5000);
methods = [true; true; true; true; false];
simulation_mc(model, methods, parameters);


%%
model = 40;
parameters = struct('sim_T',100,'MFH_w_ss',5);
methods = [true; true; true; true];
simulation_mc(model,methods,parameters)


%% PMHE1 comparison - model 42
model = 42;
parameters = struct('sim_T',100,'MFH_w_ss',2);
methods = [true; true; false; true];
simulation_mc(model,methods,parameters)

%% PMHE1 comparison - model 42
model = 42;
parameters = struct('sim_T',200,'MFH_w_ss',2,'PMHE1_N',2);
methods = [true; true; false; true; true];
simulation(model, methods, parameters);

model = 42;
parameters = struct('sim_T',200,'MFH_w_ss',2,'PMHE1_N',5);
methods = [true; true; false; true; true];
simulation(model, methods, parameters);

model = 41;
parameters = struct('sim_T',200,'PMHE1_N',5);
methods = [true; true; false; false; true];
simulation(model, methods, parameters);



%% Plots 

simN2 = load('./models/model42/simulation_N2.mat');
simN5 = load('./models/model42/simulation_N5.mat');

error_norm = zeros(5,201);
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
p = plot(0:200,error_norm,'LineWidth',3);
ylabel("$||\hat{\mathbf{x}}(k|k|k)-\mathbf{x}(k)||_2$",'Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend({'C','OS','MFH ($W_{ss} = 2$)', 'PMHE1 ($W_{\mathrm{PMHE1}} = 2$)', 'PMHE1 ($W_{\mathrm{PMHE1}} = 5$)'},'Interpreter','latex');
hold off;

mean(error_norm')




