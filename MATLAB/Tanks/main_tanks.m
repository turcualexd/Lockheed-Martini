clear, clc, close all;

B = 2;
OF = 2.42;
k_He = 5/3;

rho_f = 807;                            % kg/m^3 (289 K,  Sutton, pg. 246)
rho_ox = 1140;                          % kg/m^3 (90.4 K, Sutton, pg. 246)
m_dot_prop = 0.284;                     % kg/s

D = 1;                                  % m
H = 2;                                  % m
L_c = 0.15;                             % m
V_max = 0.8 * D^2/4 * pi * (H - L_c);   % m^3

V_He_in = V_max / B^(1/k_He);
V_prop_in = V_max - V_He_in;

M_f_in    = (1/rho_f + OF/rho_ox)^-1 * V_prop_in;
M_ox_in   = M_f_in * OF;
M_prop_in = M_f_in + M_ox_in;

V_f_in = M_f_in / rho_f;
V_ox_in = M_ox_in / rho_ox;

tb = M_prop_in / m_dot_prop / 60