clear; close all; clc;

load("cooling_cstar_cost.mat");

figure;
semilogy(axial*1e2, q_flow_temp(:,1)./1e3,"LineWidth",1.5);
hold on;
semilogy(axial*1e2, q_flow_temp(:,30)./1e3,"LineWidth",1.5);
semilogy(axial*1e2, q_flow_temp(:,end)./1e3,"LineWidth",1.5);
xlabel("$Axial\;distance\;nozzle\,[cm]$","Interpreter","latex","FontSize",15); ylabel("$\dot{q}\;[kW/m^2]$","Interpreter","latex","FontSize",15);
grid on;
xline(L_con*1e2, "LineStyle","--","LineWidth",1.5);
legend("Initial time - 0 s", "Intermediate time - 1740 s", "Final time - 3660 s","Throat section","Interpreter","latex","FontSize",14);
title("Nozzle heat flux at different time instants","Interpreter","latex","FontSize",15);
dT(1)
dT(30)
dT(62)
%%
clear; close all; clc;
load("deltaT_2.mat");

figure;
axial = axial_distance + L_con;
axial = [linspace(0,L_con,100)'; axial];
semilogy(axial*1e2, q_flow_temp(:,1)./1e3,"LineWidth",1.5);
hold on;
semilogy(axial*1e2, q_flow_temp(:,30)./1e3,"LineWidth",1.5);
semilogy(axial*1e2, q_flow_temp(:,end)./1e3,"LineWidth",1.5);
xlabel("$Axial\;distance\;nozzle\,[cm]$","Interpreter","latex","FontSize",15); ylabel("$\dot{q}\;[kW/m^2]$","Interpreter","latex","FontSize",15);
grid on;
xline(L_con*1e2, "LineStyle","--","LineWidth",1.5);
legend("Initial time - 0 s", "Intermediate time - 1740 s", "Final time - 3660 s","Throat section","Interpreter","latex","FontSize",14);
title("Nozzle heat flux at different time instants","Interpreter","latex","FontSize",15);
dT(1)
dT(30)
dT(62)