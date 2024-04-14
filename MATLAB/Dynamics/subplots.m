clear; close all; clc;
%carico file
load("mat_inc.mat");

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
subplot(3, 3, 1)
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