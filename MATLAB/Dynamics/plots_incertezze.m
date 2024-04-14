clear; close all; clc;
%carico file
load("mat_inc.mat");
 creo tvet
dt   = 10;
t_max = 4e3;
tvet = 0 : dt : t_max;
n_simulations = size(c_star_mat,1);

%calcolo vettori media

m_f_avg = mean(m_f_mat);
m_ox_avg = mean(m_ox_mat);
u_feed_f_avg = mean(u_feed_f_mat);
u_feed_ox_avg = mean(u_feed_ox_mat);
V_p_f_avg = mean(V_p_f_mat);
V_p_ox_avg =mean(V_p_ox_mat);
p_f_avg = mean(p_f_mat);
p_ox_avg = mean(p_ox_mat);
p_c_avg = mean(p_c_mat);
OF_avg = mean(OF_mat);
T_c_avg = mean(T_c_mat);
gamma_avg = mean(gamma_mat);
c_t_avg = mean(c_t_mat);
c_star_avg = mean(c_star_mat);
T_avg = mean(T_mat);
I_sp_avg = mean(I_sp_mat);
T_f_avg = mean(T_f_mat);
T_ox_avg = mean(T_ox_mat);

 
tvet=tvet-dt;
grayColor = [.6 .6 .6];
fs_t = 12;      % Titolo
fs_ax = 14;     % Assi
fs_leg = 12;    % Legenda
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),OF_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), OF_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),OF_avg(2:end), 'LineWidth',4, 'Color','r')
title("O/F ratio", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("O/F  [-]", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
%%

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),p_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), p_ox_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),p_ox_avg(2:end), 'LineWidth',1, 'Color','r')
title("Oxidizer Tank Pressure", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$P_{ox} \; [Pa]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
axes('position',[.32 .63 .25 .25])
box on
hold on
grid minor
indexofinterest= (tvet>50) & (tvet<80);
for q = 1:n_simulations
    plot(tvet(indexofinterest), p_ox_mat(q,indexofinterest), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(indexofinterest),p_ox_avg(indexofinterest), 'LineWidth',4, 'Color','r')
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),p_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), p_f_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),p_f_avg(2:end), 'LineWidth',4, 'Color','r')
% title("Fuel pressure")
% xlabel("t [s]")
% ylabel("p_{f} [Pa]")
% legend('avg', 'sim')
title("Fuel Tank Pressure", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$P_{f} \; [Pa]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
axes('position',[.28 .63 .25 .25])
box on
hold on
grid minor
indexofinterest= (tvet>50) & (tvet<80);
for q = 1:n_simulations
    plot(tvet(indexofinterest), p_f_mat(q,indexofinterest), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(indexofinterest),p_f_avg(indexofinterest), 'LineWidth',4, 'Color','r')
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),p_c_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), p_c_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),p_c_avg(2:end), 'LineWidth',4, 'Color','r')
title("Chamber Pressure", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$P_{c} \; [Pa]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
axes('position',[.27 .63 .25 .25])
box on
hold on
grid minor
indexofinterest= (tvet>50) & (tvet<100);
for q = 1:n_simulations
    plot(tvet(indexofinterest), p_c_mat(q,indexofinterest), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(indexofinterest),p_c_avg(indexofinterest), 'LineWidth',4, 'Color','r')
%%

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),m_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), m_ox_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),m_ox_avg(2:end), 'LineWidth',4, 'Color','r')
title("Oxidizer mass flow rate", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$\dot{m}_{ox} \; [kg/s]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
axes('position',[.27 .65 .25 .25])
box on
hold on
grid minor
indexofinterest= (tvet>50) & (tvet<100);
for q = 1:n_simulations
    plot(tvet(indexofinterest), m_ox_mat(q,indexofinterest), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(indexofinterest),m_ox_avg(indexofinterest), 'LineWidth',4, 'Color','r')
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),m_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), m_f_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),m_f_avg(2:end), 'LineWidth',4, 'Color','r')
title("Fuel mass flow rate", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$\dot{m}_{f} \; [kg/s]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
axes('position',[.27 .65 .25 .25])
box on
hold on
grid minor
indexofinterest= (tvet>50) & (tvet<100);
for q = 1:n_simulations
    plot(tvet(indexofinterest), m_f_mat(q,indexofinterest), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(indexofinterest),m_f_avg(indexofinterest), 'LineWidth',4, 'Color','r')
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),u_feed_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), u_feed_ox_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),u_feed_ox_avg(2:end), 'LineWidth',4, 'Color','r')
title("Oxidizer feed velocity", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$u_{feed,ox} [m/s]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),u_feed_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), u_feed_f_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
 plot(tvet(2:end),u_feed_f_avg(2:end), 'LineWidth',4, 'Color','r')
% title("Fuel feed velocity")
% xlabel("t [s]")
% ylabel("u_{feed,f} [m/s]")
% legend('avg', 'sim')
title("Fuel feed velocity", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$u_{feed,f} [m/s]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet,I_sp_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), I_sp_mat(q, (2:end)), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet,I_sp_avg, 'LineWidth',4, 'Color','r')
title("Specific Impulse", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$I_{sp} [s]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),T_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), T_ox_mat(q, (2:end)), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),T_ox_avg(2:end), 'LineWidth',4, 'Color','r')
title("Oxidizer temperature", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$T_{ox} [K]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),T_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), T_f_mat(q, (2:end)), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),T_f_avg(2:end), 'LineWidth',4, 'Color','r')
title("Fuel temperature", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$T_{f} [K]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),T_c_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), T_c_mat(q, 2:end), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),T_c_avg(2:end), 'LineWidth',4, 'Color','r')
title("Combustion Chamber temperature", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$T_{c} [K]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
%%
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),T_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), T_mat(q, 2:end), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),T_avg(2:end), 'LineWidth',2, 'Color','r')
title("Thrust", "Interpreter", "latex", "FontSize", fs_t)
xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$\mathcal{T} \; [N]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Average Value", "Simulations", "Interpreter", "latex", "FontSize", fs_leg)
axes('position',[.25 .65 .25 .25])
box on
hold on
grid minor
indexofinterest= (tvet>50) & (tvet<100);
for q = 1:n_simulations
    plot(tvet(indexofinterest), T_mat(q,indexofinterest), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(indexofinterest),T_avg(indexofinterest), 'LineWidth',4, 'Color','r')
%%