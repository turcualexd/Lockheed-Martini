clear; close all; clc; 

%dati iniziali noti

T_i_n = 1000;       % N (Aggiungere perdite per ugello)
p_c_i = 50e5;       % Pa
p_c_min = 20e5;     % Pa
rho_f = 807;        % kg/m^3
rho_ox = 1140;      % kg/m^3
L_star = 1.143;     % m
h = 2;              % m
d = 1;              % m
g0 = 9.81;          % m/s^2
R = 8.314;          % J/molK
M_m_He = 4e-3;      % kg/mol
M_m_N = 28e-3;      % kg/mol
t_max = 10000;      % s
mu_f = 0.75e-3;     % Pas
mu_ox = 0.196e-3;   % Pas

% Dati assunti
OF_i  = 2.56;        % -
eps   = 300;          % -
eps_c = 10;         % -
A_max = 0.25*pi;    % m^2
C_d = 0.7;          % -
alpha = 0.2;        % -
d_inj_f = 1e-3;     % m
d_feed = 5e-3;      % m
f_f = 0.025;         % -
f_ox = 0.015;        % -
dt = 1;             % s
lambda = 0.9986;    % -
k_He = 5/3;         % -
k_N = 7/5;          % -
T_N_f_i = 290;      % K
T_He_ox_i = 90;     % K

output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_i,'sup',eps,'case','Cazzo','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
%%
T_i = T_i_n/lambda;
V_tot = pi*d^2*h/4;
V_u = 0.8*V_tot;
T_c_i = output.froz.temperature(1);
c_star_i = output.froz.cstar(end);
c_t_i = output.froz.cf(end);
I_sp_i = c_t_i*c_star_i/g0;
m_p = T_i/(I_sp_i*g0);
A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));




