clear, clc, close all

load("final_run.mat")

% Font size
fs_t = 12;      % Titolo
fs_ax = 14;     % Assi
fs_leg = 12;    % Legenda
lw = 1;         % Linewidth

% Volumi pressurizzanti
if 0
    figure
    plot(tvet, V_p_f, 'r', LineWidth=lw)
    hold on
    grid minor
    plot(tvet, V_p_ox, 'b', LineWidth=lw)
    title("Pressurizer volumes", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$V_{pr} \; [m^3]$", "Interpreter", "latex", "FontSize", fs_ax)
    legend("Fuel branch", "Oxidizer branch", "Interpreter", "latex", "FontSize", fs_leg)
    ylim([0.1 0.9])
end

% Portate
if 0
    figure
    plot(tvet, m_f, 'r', LineWidth=lw)
    hold on
    grid minor
    plot(tvet, m_ox, 'b', LineWidth=lw)
    title("Propellant flowrates", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$\dot{m} \; [kg/s]$", "Interpreter", "latex", "FontSize", fs_ax)
    legend("Fuel", "Oxidizer", "Interpreter", "latex", "FontSize", fs_leg)
    ylim([0.02 0.22])
end

% OF
if 0
    figure
    plot(tvet, OF, 'r', LineWidth=lw)
    grid minor
    title("O/F Ratio", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("O/F $[-]$", "Interpreter", "latex", "FontSize", fs_ax)
end

% c_t
if 0
    figure
    plot(tvet, c_t, 'r', LineWidth=lw)
    grid minor
    title("Thrust coefficient", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$c_T \; [-]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylim([1.884 1.9])
end

% c_star
if 0
    figure
    plot(tvet, c_star, 'r', LineWidth=lw)
    grid minor
    title("Characteristic velocity", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$c^* \; [m/s]$", "Interpreter", "latex", "FontSize", fs_ax)
end

% Spinta
if 1
    figure
    plot(tvet, T, 'r', LineWidth=lw)
    grid minor
    title("Thrust", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$\mathcal{T} \; [N]$", "Interpreter", "latex", "FontSize", fs_ax)
end

% Impulso specifico
if 0
    figure
    plot(tvet, I_sp, 'r', LineWidth=lw)
    grid minor
    title("Specific impulse", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$I_{sp} \; [s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylim([354.5 360])
end

% Velocità di feed
if 0
    figure
    plot(tvet, u_feed_f, 'r', LineWidth=lw)
    hold on
    grid minor
    plot(tvet, u_feed_ox, 'b', LineWidth=lw)
    title("Propellant feed velocity", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$v_{fd} \; [m/s]$", "Interpreter", "latex", "FontSize", fs_ax)
    legend("Fuel", "Oxidizer", "Interpreter", "latex", "FontSize", fs_leg)
end

% Temperatura pressurizzanti (non carino) 
if 0
    figure
    plot(tvet, T_f, 'r', LineWidth=lw)
    hold on
    grid minor
    plot(tvet, T_ox, 'b', LineWidth=lw) 
    title("Pressurizer temperatures", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time [s]", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$T_{pr} \; [K]$", "Interpreter", "latex", "FontSize", fs_ax)
    legend("Fuel branch", "Oxidizer branch", "Interpreter", "latex", "FontSize", fs_leg)
end

% Temperatura di combustione
if 0
    figure
    plot(tvet, T_c, 'r', LineWidth=lw)
    grid minor
    title("Combustion chamber temperature", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$T_c \; [s]$", "Interpreter", "latex", "FontSize", fs_ax)
end

% Pressione di combustione
if 0
    figure
    plot(tvet, p_c./1e5, 'r', LineWidth=lw)
    hold on
    grid minor
    yline(20, "k--", "LineWidth", lw)
    title("Combustion chamber pressure", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$p_c \; [bar]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylim([15 50])
end

% Pressione nei tank
if 0
    figure
    plot(tvet, p_f/1e5, 'r', LineWidth=lw)
    hold on
    grid minor
    plot(tvet, p_ox/1e5, 'b', LineWidth=lw)
    title("Propellant pressures", "Interpreter", "latex", "FontSize", fs_t)
    xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
    ylabel("$p \; [bar]$", "Interpreter", "latex", "FontSize", fs_ax)
    legend("Fuel", "Oxidizer", "Interpreter", "latex", "FontSize", fs_leg)
end