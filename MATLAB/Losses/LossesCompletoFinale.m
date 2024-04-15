clear, clc, close all

%% Input data

% Initial data
rho_f = 807;        % kg/m^3
rho_ox = 1140;      % kg/m^3
R = 8.314;          % J/molK
M_m_p_ox = 4e-3;    % kg/mol
M_m_p_f = 28e-3;    % kg/mol
mu_f = 0.75e-3;     % Pas
mu_ox = 0.196e-3;   % Pas
t_max = 4000;       % s


% Constraints
T_i_n = 1000;       % N
p_c_i = 50e5;       % Pa
p_c_min = 20e5;     % Pa
h = 2;              % m
d = 1;              % m
V_tot = pi*d^2*h/4; % m^3
% 80% of V_tot usable

% Assumptions
OF_i = 2.33;        % -
eps = 300;          % -
eps_c = 10;         % -
C_d = 0.82;         % -
L_star = 1.143;     % m
alpha = 0.2;        % -
d_feed_f = 5e-3;    % m
d_feed_ox = 7e-3;   % m
dt = 1;             % s
lambda = 1;         % -
k_ox = 5/3;         % -
k_f = 7/5;          % -
T_f_i = 300;        % K
T_ox_i = 90;        % K
B = 2.78;           % -
alpha_con = 30;     % deg
d_inj_f = 1e-3;     % m


%% Nominal sizing

T_i = T_i_n/lambda;

% CEA sizing of nozzle and CC
output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_i,'sup',eps,'case','Porco Dio','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
c_star_i = output.froz.cstar(end);
c_t_i = output.froz.cf_vac(end);
T_c_i = output.froz.temperature(1);
gamma_i = output.froz.gamma(1);
I_sp_i = output.froz.isp_vac(end);
m_p_i = T_i/(c_t_i*c_star_i);

A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
A_e = eps*A_t;
A_c = eps_c*A_t;

d_e = 2*sqrt(A_e/pi);
d_t = 2*sqrt(A_t/pi);
d_c = 2*sqrt(A_c/pi);

L_con = (d_c - d_t)/(2*tand(alpha_con));

L_c = L_star/eps_c;

% Tanks sizing
V_loss = 0.25*pi*((L_c + L_con)*d^2 - L_c*d_c^2 - L_con*(d_c^2 + d_c*d_t + d_t^2)/3);

if V_loss/V_tot > 0.2
    h_t = h - L_c - L_con;
else
    h_t = h - L_c - L_con - 4*(0.2*V_tot - V_loss)/(pi*d^2);
end

V_tank_tot = pi*d^2*h_t/4;

M_f = V_tank_tot/(OF_i*(1 + 1/(B^(1/k_ox) - 1))/rho_ox + (1 + 1/(B^(1/k_f) - 1))/rho_f);
M_ox = OF_i*M_f;

V_f_i = M_f/rho_f;
V_ox_i = M_ox/rho_ox;

V_t_f = M_f*(1 + 1/(B^(1/k_f) - 1))/rho_f;
h_t_f = 4*V_t_f/(pi*d^2);

V_p_f_i = M_f/(rho_f*(B^(1/k_f) - 1));
V_p_ox_i = M_ox/(rho_ox*(B^(1/k_ox) - 1));

% Injection plate
m_f_i = m_p_i/(1 + OF_i);
m_ox_i = m_p_i*OF_i/(1 + OF_i);

dp_inj = alpha*p_c_i;

A_inj_f_tot = m_f_i/(C_d*sqrt(2*dp_inj*rho_f));
A_inj_ox_tot = m_ox_i/(C_d*sqrt(2*dp_inj*rho_ox));

A_inj_f = pi*d_inj_f^2/4;
N_inj_f = A_inj_f_tot/A_inj_f;
N_inj_f = floor(N_inj_f);

N_inj_ox = 2*N_inj_f;

A_inj_f = A_inj_f_tot/N_inj_f;
A_inj_ox = A_inj_ox_tot/N_inj_ox;

d_inj_f = 2*sqrt(A_inj_f/pi);
d_inj_ox = 2*sqrt(A_inj_ox/pi);

v_inj_f = m_f_i/(rho_f*A_inj_f_tot);
v_inj_ox = m_ox_i/(rho_ox*A_inj_ox_tot);

% Feeding lines
A_feed_f = pi*d_feed_f^2/4;
A_feed_ox = pi*d_feed_ox^2/4;

L_feed_f = h - L_c - h_t_f;
L_feed_ox = h - L_c - h_t;

u_feed_f_i = m_f_i/(rho_f*A_feed_f);
u_feed_ox_i = m_ox_i/(rho_ox*A_feed_ox);

Re_f = rho_f*u_feed_f_i*d_feed_f/mu_f;
Re_ox = rho_ox*u_feed_ox_i*d_feed_ox/mu_ox;

f_f = moody(Re_f);
f_ox = moody(Re_ox);

% Pressure losses cascade
K_f =  1 + f_f*L_feed_f/d_feed_f + (A_feed_f/(A_inj_f_tot*C_d))^2;
K_ox = 1 + f_ox*L_feed_ox/d_feed_ox + (A_feed_ox/(A_inj_ox_tot*C_d))^2;

dp_f = 0.5*rho_f*u_feed_f_i^2*K_f;
dp_ox = 0.5*rho_ox*u_feed_ox_i^2*K_ox;

p_f_i = p_c_i + dp_f;
p_ox_i = p_c_i + dp_ox;


%% Dynamics

tvet = 0 : dt : t_max;
m_f = [m_f_i nan(1, length(tvet) - 1)];
m_ox = [m_ox_i nan(1, length(tvet) - 1)];
u_feed_f = [u_feed_f_i nan(1, length(tvet) - 1)];
u_feed_ox = [u_feed_ox_i nan(1, length(tvet) - 1)];
V_p_f = [V_p_f_i nan(1, length(tvet) - 1)];
V_p_ox = [V_p_ox_i nan(1, length(tvet) - 1)];
p_f = [p_f_i nan(1, length(tvet) - 1)];
p_ox = [p_ox_i nan(1, length(tvet) - 1)];
p_c = [p_c_i nan(1, length(tvet) - 1)];
OF = [OF_i nan(1, length(tvet) - 1)];
T_c = [T_c_i nan(1, length(tvet) - 1)];
gamma = [gamma_i nan(1, length(tvet) - 1)];
c_t = [c_t_i nan(1, length(tvet) - 1)];
c_star = [c_star_i nan(1, length(tvet) - 1)];
T_id = [T_i*lambda nan(1, length(tvet) - 1)];
I_sp_id = [I_sp_i nan(1, length(tvet) - 1)];
T_f = [T_f_i nan(1, length(tvet) - 1)];
T_ox = [T_ox_i nan(1, length(tvet) - 1)];
m_p_id=[m_f_i+m_ox_i nan(1, length(tvet) - 1)];

V_f = 0;
V_ox = 0;
cont = 1;
j = 1; % p_c_it inferiore di bisezione
toll = 0.1;

while true
    
    dV_f = m_f(cont)*dt/rho_f;
    dV_ox = m_ox(cont)*dt/rho_ox;

    V_p_f_new = V_p_f(cont) + dV_f;
    V_p_ox_new = V_p_ox(cont) + dV_ox;
    
    T_f_new = T_f(cont)*(V_p_f(cont)/V_p_f_new)^(k_f - 1);
    T_ox_new = T_ox(cont)*(V_p_ox(cont)/V_p_ox_new)^(k_ox - 1);

    p_f_new = p_f(cont)*(V_p_f(cont)/V_p_f_new)^k_f;
    p_ox_new = p_ox(cont)*(V_p_ox(cont)/V_p_ox_new)^k_ox;

    V_f = V_f + dV_f;
    V_ox = V_ox + dV_ox;
    if V_f > V_f_i || V_ox > V_ox_i
        disp("Termine per volume occupato massimo raggiunto")
        break
    end
    
    p_c_new = p_c(cont);
    u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
    u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
    
    m_f_new = rho_f*A_feed_f*u_feed_f_new;
    m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
    OF_new = m_ox_new/m_f_new;
    
    output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
    
    c_star_cea = output.froz.cstar(1);
    c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
    
    while c_star_new > c_star_cea 
        
        j = j - 0.01;
        p_c_new = j*p_c(cont);

        u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
        u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
        
        m_f_new = rho_f*A_feed_f*u_feed_f_new;
        m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
        OF_new = m_ox_new/m_f_new;
        
        output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        c_star_cea = output.froz.cstar(1);
        c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
        
    end
    
    p_c_up = (j + 0.01)*p_c(cont);
    p_c_dw = j*p_c(cont);
    j = 1;

    while true

        p_c_new = (p_c_up + p_c_dw)/2;

        u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
        u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
        
        m_f_new = rho_f*A_feed_f*u_feed_f_new;
        m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
        OF_new = m_ox_new/m_f_new;
        
        output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        
        c_star_cea = output.froz.cstar(1);
        c_t_new = output.froz.cf_vac(end);
        T_c_new = output.froz.temperature(1);
        gamma_new = output.froz.gamma(1);
        c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
        T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star_new;
        I_sp_new = output.froz.isp_vac(end);
        err = abs(c_star_cea - c_star_new);
        
        if abs(err) < toll
            break
        elseif  c_star_new < c_star_cea
            p_c_dw = p_c_new;
        else
            p_c_up = p_c_new;
        end
    end

    if p_c_new < p_c_min
        disp("Terminato per pressione in camera troppo bassa")
        break
    end
    
    V_p_f(cont + 1) = V_p_f_new;
    V_p_ox(cont + 1) = V_p_ox_new;
    p_f(cont + 1) = p_f_new;
    p_ox(cont + 1) = p_ox_new;
    m_f(cont + 1) = m_f_new;
    m_ox(cont + 1) = m_ox_new;
    OF(cont + 1) = OF_new;
    T_c(cont + 1) = T_c_new;
    p_c(cont + 1) = p_c_new;
    c_t(cont + 1) = c_t_new;
    c_star(cont + 1) = c_star_cea;
    T_id(cont + 1) = T_new;
    I_sp_id(cont + 1) = I_sp_new;
    u_feed_f(cont + 1) = u_feed_f_new;
    u_feed_ox(cont + 1) = u_feed_ox_new;
    T_f(cont + 1) = T_f_new;
    T_ox(cont + 1) = T_ox_new;
    gamma(cont + 1) = gamma_new;
    m_p_id(cont+1)=m_ox_new+m_f_new;

    cont = cont + 1
end

% Rimuovi NaN
V_p_f = V_p_f(~isnan(V_p_f));
V_p_ox = V_p_ox(~isnan(V_p_ox));
p_f = p_f(~isnan(V_p_ox));
p_ox = p_ox(~isnan(V_p_ox));
m_f = m_f(~isnan(V_p_ox));
m_ox = m_ox(~isnan(V_p_ox));
OF = OF(~isnan(V_p_ox));
T_c = T_c(~isnan(V_p_ox));
p_c = p_c(~isnan(V_p_ox));
c_t = c_t(~isnan(V_p_ox));
c_star = c_star(~isnan(V_p_ox));
T_id = T_id(~isnan(V_p_ox));
I_sp_id = I_sp_id(~isnan(V_p_ox));
u_feed_f = u_feed_f(~isnan(V_p_ox));
u_feed_ox = u_feed_ox(~isnan(V_p_ox));
T_f = T_f(~isnan(V_p_ox));
T_ox = T_ox(~isnan(V_p_ox));
gamma = gamma(~isnan(V_p_ox));
tvet = tvet(1:length(gamma));

A_t_iniziale=ones(tvet);
A_t_iniziale(:)=A_t;
p_c_id=p_c;
%% Interpretation of data

dp_c_end = p_c_new - p_c_min;
I_tot_id = sum(T_id)*dt;
V_occ = (V_f - dV_f + V_ox - dV_ox + V_p_f_i + V_p_ox_i)/V_tank_tot*100;

R_p_f = R/M_m_p_f;
R_p_ox = R/M_m_p_ox;

M_p_f = p_f_i*V_p_f_i/(R_p_f*T_f_i);
M_p_ox = p_ox_i*V_p_ox_i/(R_p_ox*T_ox_i);

%% Input data

% Initial data
rho_f = 807;        % kg/m^3
rho_ox = 1140;      % kg/m^3
R = 8.314;          % J/molK
M_m_p_ox = 4e-3;    % kg/mol
M_m_p_f = 28e-3;    % kg/mol
mu_f = 0.75e-3;     % Pas
mu_ox = 0.196e-3;   % Pas
t_max = 4000;       % s


% Constraints
T_i_n = 1000;       % N
p_c_i = 50e5;       % Pa
p_c_min = 20e5;     % Pa
h = 2;              % m
d = 1;              % m
V_tot = pi*d^2*h/4; % m^3
% 80% of V_tot usable

% Assumptions
OF_i = 2.33;        % -
eps = 300;          % -
eps_c = 10;         % -
C_d = 0.82;         % -
L_star = 1.143;     % m
alpha = 0.2;        % -
d_feed_f = 5e-3;    % m
d_feed_ox = 7e-3;   % m
dt = 1;             % s
lambda = 1;         % -
k_ox = 5/3;         % -
k_f = 7/5;          % -
T_f_i = 300;        % K
T_ox_i = 90;        % K
B = 2.78;           % -
alpha_con = 30;     % deg
d_inj_f = 1e-3;     % m


%% Nominal sizing

T_i = T_i_n/lambda;

% CEA sizing of nozzle and CC
output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_i,'sup',eps,'case','Porco Dio','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
c_star_i = output.froz.cstar(end);
c_t_i = output.froz.cf_vac(end);
T_c_i = output.froz.temperature(1);
gamma_i = output.froz.gamma(1);
I_sp_i = output.froz.isp_vac(end);
m_p_i = T_i/(c_t_i*c_star_i);

A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
A_e = eps*A_t;
A_c = eps_c*A_t;

d_e = 2*sqrt(A_e/pi);
d_t = 2*sqrt(A_t/pi);
d_c = 2*sqrt(A_c/pi);

L_con = (d_c - d_t)/(2*tand(alpha_con));

L_c = L_star/eps_c;

%Initial Nozzle losses
r_t_i=sqrt(A_t/pi); % m
gamma_t_i=output.froz.gamma(2);
rho_t_i=output.froz.density(2);
mhu_t_i=output.froz.viscosity(2)*10^-4; %viscosità  
mm_t=output.froz.mw(2); %massa molora propelente in gola 
a_t_i=sqrt(gamma_t_i*output.froz.temperature(2)*(8.314/mm_t)); %velocità del suono alla throat  
r_c_i=0.382*r_t_i;   %curvatura alla gola  
Re_t_i= rho_t_i*2*r_t_i*a_t_i/mhu_t_i; %reynolds in gola 
Re_mod_i=sqrt((r_t_i)/r_c_i)*Re_t_i; %Reynolds modificato 
%discharge coefficient 
Cd_i=1-(((gamma_t_i+1)/2)^(3/4))*(3.266-(2.128/(gamma_t_i+1)))*(Re_mod_i^(-1/2))+0.9428*((gamma_t_i-1)*(gamma_t_i+2)/((gamma_t_i+1)^(1/2)))*(Re_mod_i^-1); 
%cd molto un po' più alto dei valori tipici dati da maggi  
m_re_i=Cd_i*m_p_i;
% Area di gola relativa 
A_t_relativa_i=(c_star_i*m_re_i)/(output.froz.pressure(1)*10^5); %stare attendi che output di Cea è in bar 

%T e Is iiniziale reali
T_re_i=lambda*(m_re_i)*c_t_i*c_star_i;
I_sp_re_i= T_re_i/(m_re_i*9.81);

% Tanks sizing
V_loss = 0.25*pi*((L_c + L_con)*d^2 - L_c*d_c^2 - L_con*(d_c^2 + d_c*d_t + d_t^2)/3);

if V_loss/V_tot > 0.2
    h_t = h - L_c - L_con;
else
    h_t = h - L_c - L_con - 4*(0.2*V_tot - V_loss)/(pi*d^2);
end

V_tank_tot = pi*d^2*h_t/4;

M_f = V_tank_tot/(OF_i*(1 + 1/(B^(1/k_ox) - 1))/rho_ox + (1 + 1/(B^(1/k_f) - 1))/rho_f);
M_ox = OF_i*M_f;

V_f_i = M_f/rho_f;
V_ox_i = M_ox/rho_ox;

V_t_f = M_f*(1 + 1/(B^(1/k_f) - 1))/rho_f;
h_t_f = 4*V_t_f/(pi*d^2);

V_p_f_i = M_f/(rho_f*(B^(1/k_f) - 1));
V_p_ox_i = M_ox/(rho_ox*(B^(1/k_ox) - 1));

% Injection plate
m_f_i = m_p_i/(1 + OF_i);
m_ox_i = m_p_i*OF_i/(1 + OF_i);

dp_inj = alpha*p_c_i;

A_inj_f_tot = m_f_i/(C_d*sqrt(2*dp_inj*rho_f));
A_inj_ox_tot = m_ox_i/(C_d*sqrt(2*dp_inj*rho_ox));

A_inj_f = pi*d_inj_f^2/4;
N_inj_f = A_inj_f_tot/A_inj_f;
N_inj_f = floor(N_inj_f);

N_inj_ox = 2*N_inj_f;

A_inj_f = A_inj_f_tot/N_inj_f;
A_inj_ox = A_inj_ox_tot/N_inj_ox;

d_inj_f = 2*sqrt(A_inj_f/pi);
d_inj_ox = 2*sqrt(A_inj_ox/pi);

v_inj_f = m_f_i/(rho_f*A_inj_f_tot);
v_inj_ox = m_ox_i/(rho_ox*A_inj_ox_tot);

% Feeding lines
A_feed_f = pi*d_feed_f^2/4;
A_feed_ox = pi*d_feed_ox^2/4;

L_feed_f = h - L_c - h_t_f;
L_feed_ox = h - L_c - h_t;

u_feed_f_i = m_f_i/(rho_f*A_feed_f);
u_feed_ox_i = m_ox_i/(rho_ox*A_feed_ox);

Re_f = rho_f*u_feed_f_i*d_feed_f/mu_f;
Re_ox = rho_ox*u_feed_ox_i*d_feed_ox/mu_ox;

f_f = moody(Re_f);
f_ox = moody(Re_ox);

% Pressure losses cascade
K_f =  1 + f_f*L_feed_f/d_feed_f + (A_feed_f/(A_inj_f_tot*C_d))^2;
K_ox = 1 + f_ox*L_feed_ox/d_feed_ox + (A_feed_ox/(A_inj_ox_tot*C_d))^2;

dp_f = 0.5*rho_f*u_feed_f_i^2*K_f;
dp_ox = 0.5*rho_ox*u_feed_ox_i^2*K_ox;

p_f_i = p_c_i + dp_f;
p_ox_i = p_c_i + dp_ox;

% modifiche per m_ox e m_fu
m_ox_new=(OF_i/(1+OF_i))*m_re_i;
m_fu_new=m_re_i-m_ox_i;
%% Dynamics

tvet = 0 : dt : t_max;
m_f = [m_f_i nan(1, length(tvet) - 1)];
m_ox = [m_ox_i nan(1, length(tvet) - 1)];
u_feed_f = [u_feed_f_i nan(1, length(tvet) - 1)];
u_feed_ox = [u_feed_ox_i nan(1, length(tvet) - 1)];
V_p_f = [V_p_f_i nan(1, length(tvet) - 1)];
V_p_ox = [V_p_ox_i nan(1, length(tvet) - 1)];
p_f = [p_f_i nan(1, length(tvet) - 1)];
p_ox = [p_ox_i nan(1, length(tvet) - 1)];
p_c = [p_c_i nan(1, length(tvet) - 1)];
OF = [OF_i nan(1, length(tvet) - 1)];
T_c = [T_c_i nan(1, length(tvet) - 1)];
gamma = [gamma_i nan(1, length(tvet) - 1)];
c_t = [c_t_i nan(1, length(tvet) - 1)];
c_star = [c_star_i nan(1, length(tvet) - 1)];
T = [T_i*lambda nan(1, length(tvet) - 1)];
I_sp = [I_sp_i nan(1, length(tvet) - 1)];
T_f = [T_f_i nan(1, length(tvet) - 1)];
T_ox = [T_ox_i nan(1, length(tvet) - 1)];

%Nozzle losses variable initialisation
r_t = [r_t_i nan(1, length(tvet) - 1)];
gamma_t= [gamma_t_i nan(1, length(tvet) - 1)];
rho_t= [rho_t_i nan(1, length(tvet) - 1)];
mhu_t= [mhu_t_i nan(1, length(tvet) - 1)];
a_t= [a_t_i nan(1, length(tvet) - 1)];  
r_c= [r_c_i nan(1, length(tvet) - 1)]; 
Re_t= [Re_t_i nan(1, length(tvet) - 1)]; 
Re_mod= [Re_mod_i nan(1, length(tvet) - 1)]; 
Cd_Ale= [Cd_i nan(1, length(tvet) - 1)]; 
m_re= [m_re_i nan(1, length(tvet) - 1)];
E_rate=3.1886*10^-8; %erosion rate [m/s]
m_p= [m_p_i nan(1, length(tvet) - 1)];
A_t_sal= [A_t nan(1, length(tvet) - 1)];
A_t_relativa = [A_t_relativa_i nan(1, length(tvet) - 1)];
T_re = [T_re_i nan(1, length(tvet) - 1)];
I_sp_re = [I_sp_re_i nan(1, length(tvet) - 1)];

V_f = 0;
V_ox = 0;
cont = 1;
j = 1; % p_c_it inferiore di bisezione
toll = 0.1;

while true

    %nozzle erosion calculation
    r_t(cont+1)=r_t(cont)+ E_rate;
    A_t=pi*(r_t(cont+1)^2);

    dV_f = m_f(cont)*dt/rho_f;
    dV_ox = m_ox(cont)*dt/rho_ox;

    V_p_f_new = V_p_f(cont) + dV_f;
    V_p_ox_new = V_p_ox(cont) + dV_ox;
    
    T_f_new = T_f(cont)*(V_p_f(cont)/V_p_f_new)^(k_f - 1);
    T_ox_new = T_ox(cont)*(V_p_ox(cont)/V_p_ox_new)^(k_ox - 1);

    p_f_new = p_f(cont)*(V_p_f(cont)/V_p_f_new)^k_f;
    p_ox_new = p_ox(cont)*(V_p_ox(cont)/V_p_ox_new)^k_ox;

    V_f = V_f + dV_f;
    V_ox = V_ox + dV_ox;
    if V_f > V_f_i || V_ox > V_ox_i
        disp("Termine per volume occupato massimo raggiunto")
        break
    end
    
    p_c_new = p_c(cont);
    u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
    u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
    
    m_f_new = rho_f*A_feed_f*u_feed_f_new;
    m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
    OF_new = m_ox_new/m_f_new;
    
    output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
    
    c_star_cea = output.froz.cstar(1);
    c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
    
    while c_star_new > c_star_cea 
        
        j = j - 0.01;
        p_c_new = j*p_c(cont);

        u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
        u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
        
        m_f_new = rho_f*A_feed_f*u_feed_f_new;
        m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
        OF_new = m_ox_new/m_f_new;
        
        output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        c_star_cea = output.froz.cstar(1);
        c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
        
    end
    
    p_c_up = (j + 0.01)*p_c(cont);
    p_c_dw = j*p_c(cont);
    j = 1;

    while true

        p_c_new = (p_c_up + p_c_dw)/2;

        u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
        u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
        
        m_f_new = rho_f*A_feed_f*u_feed_f_new;
        m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
        OF_new = m_ox_new/m_f_new;
        
        output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        
        c_star_cea = output.froz.cstar(1);
        c_t_new = output.froz.cf_vac(end);
        T_c_new = output.froz.temperature(1);
        gamma_new = output.froz.gamma(1);
        c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
        T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star_new;
        I_sp_new = output.froz.isp_vac(end);
        err = abs(c_star_cea - c_star_new);
        
        if abs(err) < toll
            break
        elseif  c_star_new < c_star_cea
            p_c_dw = p_c_new;
        else
            p_c_up = p_c_new;
        end
    end
    
    c_star_new = c_star_cea;

    if p_c_new < p_c_min
        disp("Terminato per pressione in camera troppo bassa")
        break
    end
   
    % Computation of boundary layer losses 
    gamma_t_new=output.froz.gamma(2);
    rho_t_new=output.froz.density(2);
    mhu_t_new=output.froz.viscosity(2)*10^-4; %viscosità  
    mm_t=output.froz.mw(2); %massa molora propelente in gola 
    a_t_new=sqrt(gamma_t_new*output.froz.temperature(2)*(8.314/mm_t)); %velocità del suono alla throat  
    r_c_new=0.382*r_t(1);   %curvatura alla gola  
    Re_t_new= rho_t_new*r_t(1)*2*a_t_new/mhu_t_new; %reynolds in gola 
    Re_mod_new=sqrt((r_t(1))/r_c_new)*Re_t_new; %Reynolds modificato 
    %discharge coefficient 
    Cd_new=1-(((gamma_t_new+1)/2)^(3/4))*(3.266-(2.128/(gamma_t_new+1)))*(Re_mod_new^(-1/2))+0.9428*((gamma_t_new-1)*(gamma_t_i+2)/((gamma_t_new+1)^(1/2)))*(Re_mod_new^-1); 
    % mi trovo m_p con m fuel + m ox
    m_re_new=Cd_new*(m_f_new + m_ox_new);
    A_t_relativa_new=(c_star_cea*m_re_new)/(output.froz.pressure(1)*10^5);
    
    %Trust and Specific impulse for the real case 
    T_re_new = lambda*(m_re_new)*c_t_new*c_star_new;
    I_sp_re_new=T_re_new/(m_re_new*9.81);

    m_ox_new=(OF_new/(1+OF_new))*m_re_new;
    m_fu_new=m_re_new-m_ox_new;

    V_p_f(cont + 1) = V_p_f_new;
    V_p_ox(cont + 1) = V_p_ox_new;
    p_f(cont + 1) = p_f_new;
    p_ox(cont + 1) = p_ox_new;
    m_f(cont + 1) = m_f_new;
    m_ox(cont + 1) = m_ox_new;
    OF(cont + 1) = OF_new;
    T_c(cont + 1) = T_c_new;
    p_c(cont + 1) = p_c_new;
    c_t(cont + 1) = c_t_new;
    c_star(cont + 1) = c_star_cea;
    T(cont + 1) = T_new;
    I_sp(cont + 1) = I_sp_new;
    u_feed_f(cont + 1) = u_feed_f_new;
    u_feed_ox(cont + 1) = u_feed_ox_new;
    T_f(cont + 1) = T_f_new;
    T_ox(cont + 1) = T_ox_new;
    gamma(cont + 1) = gamma_new;

    % losses variable save
    gamma_t(cont+1)= gamma_t_new;
    rho_t(cont+1)= rho_t_new;
    mhu_t(cont+1)= mhu_t_new;
    a_t(cont+1)= a_t_new;  
    r_c(cont+1)= r_c_new; 
    Re_t(cont+1)= Re_t_new; 
    Re_mod(cont+1)= Re_mod_new; 
    Cd_Ale(cont+1)= Cd_new; 
    m_re(cont+1)= m_re_new;
    m_p(cont+1)= m_f_new + m_ox_new;
    A_t_sal(cont+1)=A_t;
    A_t_relativa(cont+1)=A_t_relativa_new;
    T_re(cont+1)=T_re_new;
    I_sp_re(cont+1)=I_sp_re_new;

    cont = cont + 1
end
%%
% Rimuovi NaN
V_p_f = V_p_f(~isnan(V_p_f));
V_p_ox = V_p_ox(~isnan(V_p_ox));
p_f = p_f(~isnan(V_p_ox));
p_ox = p_ox(~isnan(V_p_ox));
m_f = m_f(~isnan(V_p_ox));
m_ox = m_ox(~isnan(V_p_ox));
OF = OF(~isnan(V_p_ox));
T_c = T_c(~isnan(V_p_ox));
p_c = p_c(~isnan(V_p_ox));
c_t = c_t(~isnan(V_p_ox));
c_star = c_star(~isnan(V_p_ox));
T = T(~isnan(V_p_ox));
I_sp = I_sp(~isnan(V_p_ox));
u_feed_f = u_feed_f(~isnan(V_p_ox));
u_feed_ox = u_feed_ox(~isnan(V_p_ox));
T_f = T_f(~isnan(V_p_ox));
T_ox = T_ox(~isnan(V_p_ox));
gamma = gamma(~isnan(V_p_ox));
tvet = tvet(1:length(gamma));
T_re = T_re(~isnan(V_p_ox));

%% Interpretation of data

dp_c_end = p_c_new - p_c_min;
I_tot_re = sum(T_re)*dt;
V_occ = (V_f - dV_f + V_ox - dV_ox + V_p_f_i + V_p_ox_i)/V_tank_tot*100;

R_p_f = R/M_m_p_f;
R_p_ox = R/M_m_p_ox;

M_p_f = p_f_i*V_p_f_i/(R_p_f*T_f_i);
M_p_ox = p_ox_i*V_p_ox_i/(R_p_ox*T_ox_i);

%% Grafici 
clear, clc, close all;
load("LossesCompletoFinale.mat")

lw=1.5;
fs_ax = 15;
fs_t=15;
fs_leg=14;
%confronto thrust
figure 
plot(T_re, LineWidth=lw)
hold on
plot(T_id, LineWidth=lw)
title("Thrust comparison", "Interpreter", "latex", "FontSize", fs_t)
legend("Real thrust","Ideal thrust", "Interpreter", "latex", "FontSize", fs_leg)
grid minor
xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$\mathcal{T} \; [N]$", "Interpreter", "latex", "FontSize", fs_ax)

%confronto portate
figure 
plot(m_re, LineWidth=lw)
hold on
plot(m_p_id, LineWidth=lw)
title("Mass flow rate comparison", "Interpreter", "latex", "FontSize", fs_t)
legend("Real mass flow rate","Ideal mass flow rate", "Interpreter", "latex", "FontSize", fs_leg)
grid minor
xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$\dot{m} \; [kg/s]$", "Interpreter", "latex", "FontSize", fs_ax)

%confronto impulso totale
figure 
plot(I_sp_re, LineWidth=lw)
hold on 
plot(I_sp_id, LineWidth=lw)
legend("Real specific impulse","Ideal specific impulse", "Interpreter", "latex", "FontSize", fs_leg)
title("Specific impulse comparison", "Interpreter", "latex", "FontSize", fs_t)
grid minor
xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$I_{sp} \; [s]$", "Interpreter", "latex", "FontSize", fs_ax)

figure 
plot(Cd_Ale,LineWidth=lw)
legend("$C_d$", "Interpreter", "latex", "FontSize", fs_leg)
title("Discharge coefficient variation", "Interpreter", "latex", "FontSize", fs_t)
grid minor
xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
ylabel("$C_d \; [-]$", "Interpreter", "latex", "FontSize", fs_ax)


xx=1:length(T_re);
figure
grid minor
hold on
plot(A_t_sal, LineWidth=lw)
plot(A_t_relativa, LineWidth=lw)
yline(A_t_sal(1), 'k--', LineWidth=lw)
title("Throat area variation", "Interpreter", "latex", "FontSize", fs_t)
ylabel("$A \; [m^2]$", "Interpreter", "latex", "FontSize", fs_ax)
xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Real throat area", "Apparent throat area", "Nominal throat area", "Interpreter", "latex", "FontSize", fs_leg)
ylim([1e-4 1.1e-4])

figure
plot(p_c, LineWidth=lw)
hold on 
grid minor
plot(p_c_id, LineWidth=lw)
title("Chamber pressure comparison", "Interpreter", "latex", "FontSize", fs_t)
ylabel("$p_c\; [Pa]$", "Interpreter", "latex", "FontSize", fs_ax)
xlabel("time $[s]$", "Interpreter", "latex", "FontSize", fs_ax)
legend("Real chamber pressure", "Ideal chamber pressure","Interpreter", "latex", "FontSize", fs_leg)


% stampo gli impulsi totali 
I_tot_re
I_tot_id