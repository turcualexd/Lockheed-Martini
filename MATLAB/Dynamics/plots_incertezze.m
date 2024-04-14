clear; close all; clc;
%% carico file
load("mat_inc.mat");
%% creo tvet
dt   = 1;
t_max = 4e3;
tvet = 0 : dt : t_max;
%% calcolo vettori media

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

%% faccio i plots
tvet=tvet-dt;
grayColor = [.6 .6 .6];
figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),OF_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), OF_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),OF_avg(2:end), 'LineWidth',4, 'Color','r')
title("O/F Ratio")
xlabel("t [s]")
ylabel("O/F [-]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),p_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), p_ox_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),p_ox_avg(2:end), 'LineWidth',4, 'Color','r')
title("Oxidizer pressure")
xlabel("t [s]")
ylabel("p_{ox} [Pa]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),p_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), p_f_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),p_f_avg(2:end), 'LineWidth',4, 'Color','r')
title("Fuel pressure")
xlabel("t [s]")
ylabel("p_{f} [Pa]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),p_c_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), p_c_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),p_c_avg(2:end), 'LineWidth',4, 'Color','r')
title("Chamber pressure")
xlabel("t [s]")
ylabel("p_{c} [Pa]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),m_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), m_ox_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),m_ox_avg(2:end), 'LineWidth',4, 'Color','r')
title("Oxidizer mass flow rate")
xlabel("t [s]")
ylabel("m_{ox} [kg/s]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),m_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), m_f_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),m_f_avg(2:end), 'LineWidth',4, 'Color','r')
title("Fuel mass flow rate")
xlabel("t [s]")
ylabel("m_{f} [kg/s]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),u_feed_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), u_feed_ox_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),u_feed_ox_avg(2:end), 'LineWidth',4, 'Color','r')
title("Oxidizer feed velocity")
xlabel("t [s]")
ylabel("u_{feed,ox} [m/s]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),u_feed_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), u_feed_f_mat(q, (2:end)),'-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),u_feed_f_avg(2:end), 'LineWidth',4, 'Color','r')
title("Fuel feed velocity")
xlabel("t [s]")
ylabel("u_{feed,f} [m/s]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet,I_sp_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), I_sp_mat(q, (2:end)), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet,I_sp_avg, 'LineWidth',4, 'Color','r')
title("Specific Impulse")
xlabel("t [s]")
ylabel("I_{sp} [s]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),T_ox_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), T_ox_mat(q, (2:end)), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),T_ox_avg(2:end), 'LineWidth',4, 'Color','r')
title("Oxidizer temperature")
xlabel("t [s]")
ylabel("T_{oxidizer} [K]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),T_f_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), T_f_mat(q, (2:end)), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),T_f_avg(2:end), 'LineWidth',4, 'Color','r')
title("Fuel temperature")
xlabel("t [s]")
ylabel("T_{fuel} [K]")
legend('avg', 'sim')

figure
hold on
grid minor
xlim([0 3500])
plot(tvet(2:end),T_c_avg(2:end), 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet(2:end), T_c_mat(q, 2:end), '-','LineWidth',0.5, 'Color',grayColor);
end
plot(tvet(2:end),T_c_avg(2:end), 'LineWidth',4, 'Color','r')
title("Combustion Chamber temperature")
xlabel("t [s]")
ylabel("T_{cc} [K]")
legend('avg', 'sim')