%% Dimensionamento Alex
clear, clc, close all

% Dati iniziali noti
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
t_max = 100000;      % s
mu_f = 0.75e-3;     % Pas
mu_ox = 0.196e-3;   % Pas

% Dati assunti

OF_i = 2.4;        % -
% OF_final = 2.5;    % -
eps = 300;          % -
eps_c = 15;         % -
A_max = 0.25*pi;    % m^2
C_d = 0.7;          % - cambio questo
alpha = 0.05;       % - qui cambio il valore
d_inj_f = 1e-3;     % m
d_feed = 5e-3;      % m
f_f = 0.025;        % -
f_ox = 0.015;       % -
dt = 1;             % s
lambda = 0.9986;    % -
k_He = 5/3;         % -
k_N = 7/5;          % -
T_N_f_i = 300;      % K
T_He_ox_i = 90;     % K
correzione = 1;     % -
p_voluta = 35e5;    % Pa

% Dimensionamento a ritroso
T_i = T_i_n/lambda;
V_tot = pi*d^2*h/4;
V_u = 0.8*V_tot;

output = cea(CEA('problem','rkt','nfz',1,'o/f',OF_i,'sup',eps,'case','Cazzo','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
c_star_i = output.froz.cstar(end);
c_t_i = output.froz.cf(end);
T_c_i = output.froz.temperature(1);
gamma_i = output.froz.gamma(1);
I_sp_i = c_t_i*c_star_i/g0;
m_p_i = T_i/(c_t_i*c_star_i);
A_t = c_star_i * m_p_i / p_c_i;
%A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
A_e = eps*A_t;

A_c = eps_c*A_t;
d_c = 2*sqrt(A_c/pi);
L_c = L_star/eps_c;
V_int_c = pi*(d^2 - d_c^2)*L_c/4;

if V_int_c/V_tot > 0.2
    h_tank = h - L_c;
else
    h_tank = h - L_c - 4*(0.2*V_tot - V_int_c)/(pi*d^2);
end

V_tank_tot = pi*d^2*h_tank/4;
m_f_i = m_p_i/(1 + OF_i);
m_ox_i = m_p_i*OF_i/(1 + OF_i);
dp_inj = alpha*p_c_i;
A_inj_f_tot = m_f_i/(C_d*sqrt(2*dp_inj*rho_f));
A_inj_ox_tot = m_ox_i/(C_d*sqrt(2*dp_inj*rho_ox));

L = h - L_c;
L_feed_f = 2*L/3;   % Assunto
L_feed_ox = L/3;    % Assunto
A_feed = pi*d_feed^2/4;
u_feed_f_i = m_f_i/(rho_f*A_feed*correzione);
u_feed_ox_i = m_ox_i/(rho_ox*A_feed);

K_f = 1 + f_f*L_feed_f/(sqrt(correzione)*d_feed) + (A_feed*correzione/(A_inj_f_tot*C_d))^2; % qui stiamo modificando le perdite
K_ox = 1 + f_ox*L_feed_ox/d_feed + (A_feed/(A_inj_ox_tot*C_d))^2 ;

dp_f = 0.5*rho_f*u_feed_f_i^2*K_f;
dp_ox = 0.5*rho_ox*u_feed_ox_i^2*K_ox;

p_f_i = p_c_i + dp_f;
p_ox_i = p_c_i + dp_ox;

v_N_f_i = R*T_N_f_i/(M_m_N*p_f_i);
rho_N = 1/v_N_f_i;

v_He_ox_i = R*T_He_ox_i/(M_m_He*p_ox_i);
rho_He = 1/v_He_ox_i;

V_N_f_i = 0.16*V_tot
V_He_ox_i = OF_i * V_N_f_i
V_tot
V_press_i = V_N_f_i + V_He_ox_i
perce_press = (V_N_f_i + V_He_ox_i)/V_tot

% Iterazione
tvet = 0 : dt : t_max;
m_f = [m_f_i nan(1, length(tvet) - 1)];
m_ox = [m_ox_i nan(1, length(tvet) - 1)];
u_feed_f = [u_feed_f_i nan(1, length(tvet) - 1)];
u_feed_ox = [u_feed_ox_i nan(1, length(tvet) - 1)];
V_N_f = [V_N_f_i nan(1, length(tvet) - 1)];
V_He_ox = [V_He_ox_i nan(1, length(tvet) - 1)];
p_f = [p_f_i nan(1, length(tvet) - 1)];
p_ox = [p_ox_i nan(1, length(tvet) - 1)];
OF = [OF_i nan(1, length(tvet) - 1)];
T_c = [T_c_i nan(1, length(tvet) - 1)];
gamma = [gamma_i nan(1, length(tvet) - 1)];
p_c = [p_c_i nan(1, length(tvet) - 1)];
c_t = [c_t_i nan(1, length(tvet) - 1)];
c_star = [c_star_i nan(1, length(tvet) - 1)];
T = [T_i*lambda nan(1, length(tvet) - 1)];
I_sp = [I_sp_i nan(1, length(tvet) - 1)];

V_f = 0;
V_ox = 0;
valido = 1;
valido_2 = 1;
cont = 1;
j = 1; % p_c_it inferiore di bisezione
toll = 0.1;

while valido
    dV_f = m_f(cont)*dt/rho_f;
    dV_ox = m_ox(cont)*dt/rho_ox;

    V_N_f_new = V_N_f(cont) + dV_f;
    V_He_ox_new = V_He_ox(cont) + dV_ox;

    p_f_new = p_f(cont)*(V_N_f(cont)/V_N_f_new)^k_N;
    p_ox_new = p_ox(cont)*(V_He_ox(cont)/V_He_ox_new)^k_He;
    
    V_f = V_f + dV_f;
    V_ox = V_ox + dV_ox;

    if (V_f + V_ox + V_He_ox_i + V_N_f_i) > V_tank_tot
        valido = 0; % QUI DEVE ESSERCI 0
        disp("Termine per volume occupato massimo raggiunto")
        continue
    end
    
    p_c_it = p_c(cont);

    u_feed_f_it = sqrt(2*(p_f_new - p_c_it)/(rho_f*K_f));
    u_feed_ox_it = sqrt(2*(p_ox_new - p_c_it)/(rho_ox*K_ox));
    
    m_f_it = rho_f*A_feed*u_feed_f_it*correzione;
    m_ox_it = rho_ox*A_feed*u_feed_ox_it;
    OF_it = m_ox_it/m_f_it;
    
    output = cea(CEA('problem','rkt','nfz',1,'o/f',OF_it,'sup',eps,'case','Cazzo','p,bar',p_c_it/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
    c_star_cea = output.froz.cstar(1);
    c_star_it = A_t*p_c_it/(m_f_it + m_ox_it);
    
    while c_star_it > c_star_cea 
        
        j = j - 0.01;
        p_c_it = j*p_c(cont);

        u_feed_f_it = sqrt(2*(p_f_new - p_c_it)/(rho_f*K_f));
        u_feed_ox_it = sqrt(2*(p_ox_new - p_c_it)/(rho_ox*K_ox));
        
        m_f_it = rho_f*A_feed*u_feed_f_it*correzione;
        m_ox_it = rho_ox*A_feed*u_feed_ox_it;
        OF_it = m_ox_it/m_f_it;
        
        output = cea(CEA('problem','rkt','nfz',1,'o/f',OF_it,'sup',eps,'case','Cazzo','p,bar',p_c_it/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        c_star_cea = output.froz.cstar(1);
        c_star_it = A_t*p_c_it/(m_f_it + m_ox_it);
        
    end
    
    p_c_up = (j + 0.01)*p_c(cont);
    p_c_dw = j*p_c(cont);
    j = 1;

    while valido_2

        p_c_new = (p_c_up + p_c_dw)/2;

        u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
        u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
        
        m_f_new = rho_f*A_feed*u_feed_f_new*correzione;
        m_ox_new = rho_ox*A_feed*u_feed_ox_new;
        OF_new = m_ox_new/m_f_new;
        
        output = cea(CEA('problem','rkt','nfz',1,'o/f',OF_new,'sup',eps,'case','Cazzo','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        c_star_cea = output.froz.cstar(1);
        c_t_new = output.froz.cf(end);
        T_c_new = output.froz.temperature(1);
        c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
        T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star_new;
        I_sp_new = c_t_new*c_star_new/g0;
        err = abs(c_star_cea - c_star_new);
        
        if abs(err) < toll
            valido_2 = 0;
        elseif  c_star_new < c_star_cea
            p_c_dw = p_c_new;
        else
            p_c_up = p_c_new;
        end
    end
    
    valido_2 = 1;

    if p_c_new < p_c_min
        valido = 0;
        disp("Terminato per pressione in camera troppo bassa")
        perce = (V_f + V_ox + V_He_ox_i + V_N_f_i)/(V_tank_tot) * 100
        ox_mass = V_ox * rho_ox
        f_mass = V_f * rho_f
        continue
    end
    
    V_N_f(cont + 1) = V_N_f_new;
    V_He_ox(cont + 1) = V_He_ox_new;
    p_f(cont + 1) = p_f_new;
    p_ox(cont + 1) = p_ox_new;
    m_f(cont + 1) = m_f_new;
    m_ox(cont + 1) = m_ox_new;
    OF(cont + 1) = OF_new;
    T_c(cont + 1) = T_c_new;
    p_c(cont + 1) = p_c_new;
    c_t(cont + 1) = c_t_new;
    c_star(cont + 1) = c_star_new;
    T(cont + 1) = T_new;
    I_sp(cont + 1) = I_sp_new;
    u_feed_f(cont + 1) = u_feed_f_new;
    u_feed_ox(cont + 1) = u_feed_ox_new;

    cont = cont + 1
end

%% figure
plot(OF)
grid minor
title("OF")

figure
plot(p_f, 'r')
grid minor
hold on
plot(p_ox, 'b')
plot(p_c, 'm')
title("Pressioni")
legend("p_f", "p_{ox}", "p_c")

figure
plot(m_f, 'r')
hold on
grid minor
plot(m_ox, 'b')
plot(m_ox + m_f, 'g')
title("Portate")
legend("m_f", "m_{ox}", "m_{tot}")

figure
plot(u_feed_f, 'r')
hold on
grid minor
plot(u_feed_ox, 'b')
title("Velocità feed")
legend("u_{feed,f}", "u_{feed,ox}")

figure
plot(V_N_f, 'r')
hold on
grid minor
plot(V_He_ox, 'b')
title("Volumi pressurizzant")
legend("Volume N_2", "Volume He")

figure
plot(I_sp)
grid minor
title("Impulso specifico")

figure
plot(T_c)
grid minor
title("Temperatura camera [K]");

T = T(~isnan(T));
I_tot = sum(T)*dt;