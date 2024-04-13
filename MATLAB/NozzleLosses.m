% Nozzle losses 
clc,clear, close all;

%% Visosity displayment 
% 
% Dati
% rho_t=
% mhu_t=
% r_t=   %Ricavato dalla vostra At iniziale 
% gamma=
% m_id=
% 
% r_c=1/(1.5*r_t);   %curvatura alla gola  
% Re_t= rho_t*r_t/mhu_t;
% Re_mod=sqrt((r_t)/r_c)*Re_t;
% 
% %discharge coefficient 
% Cd=1-((gamma+1)/2)^(3/4)*(3.266-(2.128/(gamma+1)))*Re_mod^(-1/2)+0.9428*((gamma-1)*(gamma+2)/((gamma+1)^(1/2)))*Re_mod^-1;
% 
% m_re=Cd*m_id;

%% Throat Erosion 
%E_rate=  [mm/s] sulla letteratura con razzi con dimensione simili alle nostre 
% c'è 0,05 mm/s (Guardare Nozzle tesi)sulle slide di maggi 0,1-0,5 mm/s 
% ma non è specificato il tipo di razzo 

% Per valutare questo condizioni ho deciso di modificare cose 2 in modo che
% ad ogni iterazione At

%% C_star variabile
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
t_max = 10000;      % s
mu_f = 0.75e-3;     % Pas
mu_ox = 0.196e-3;   % Pas

% Dati assunti
OF_i = 2.42;        % -
OF_m = 2.56;        % -
eps = 300;          % -
eps_c = 10;         % -
C_d = 0.7;          % -
alpha = 0.2;        % -
d_feed_f = 5e-3;    % m
d_feed_ox = 7e-3;   % m
dt = 1;             % s
lambda = 0.9974;    % -
k_ox = 5/3;         % -
k_f = 7/5;          % -
T_f_i = 300;        % K
T_ox_i = 90;        % K
B_f = 2.6;          % -
B_ox = 2.6;         % -


% Dimensionamento a ritroso
T_i = T_i_n/lambda;
V_tot = pi*d^2*h/4;
output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_i,'sup',eps,'case','Porco Dio','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
c_star_i = output.froz.cstar(end);
c_t_i = output.froz.cf(end);
T_c_i = output.froz.temperature(1);
gamma_i = output.froz.gamma(1);
I_sp_i = output.froz.isp(end);
m_p_i = T_i/(c_t_i*c_star_i);
A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
A_e = eps*A_t;
A_c = eps_c*A_t;
d_c = 2*sqrt(A_c/pi);
L_c = L_star/eps_c;
V_int_c = pi*(d^2 - d_c^2)*L_c/4;

% Modifice Ale Viscosity displaiment 
r_t_i=sqrt(A_t/pi); % m
gamma_t_i=output.froz.gamma(2);
rho_t_i=output.froz.density(2);
mhu_t_i=output.froz.viscosity(2)*10^-4; %viscosità  
mm_t=output.froz.mw(2); %massa molora propelente in gola 
a_t_i=sqrt(gamma_t_i*output.froz.temperature(2)*(8.314/mm_t)); %velocità del suono alla throat  
r_c_i=1.5*r_t_i;   %curvatura alla gola  
Re_t_i= rho_t_i*2*r_t_i*a_t_i/mhu_t_i; %reynolds in gola 
Re_mod_i=sqrt((r_t_i)/r_c_i)*Re_t_i; %Reynolds modificato 
%discharge coefficient 
Cd_i=1-(((gamma_t_i+1)/2)^(3/4))*(3.266-(2.128/(gamma_t_i+1)))*(Re_mod_i^(-1/2))+0.9428*((gamma_t_i-1)*(gamma_t_i+2)/((gamma_t_i+1)^(1/2)))*(Re_mod_i^-1); 
%cd molto un po' più alto dei valori tipici dati da maggi  
m_re_i=Cd_i*m_p_i;

% Area di gola relativa 
A_t_relativa_i=(c_star_i*m_re_i)/(output.froz.pressure(1)*10^5); %stare attendi che output di Cea è in bar 

%%
if V_int_c/V_tot > 0.2
    h_tank = h - L_c;
else
    h_tank = h - L_c - 4*(0.2*V_tot - V_int_c)/(pi*d^2);
end

V_tank_tot = pi*d^2*h_tank/4;

% Sistema Alex
M_f = V_tank_tot/(OF_m*(1 + 1/(B_ox^(1/k_ox) - 1))/rho_ox + (1 + 1/(B_f^(1/k_f) - 1))/rho_f);
M_ox = OF_m*M_f;

V_f_i = M_f/rho_f;
V_ox_i = M_ox/rho_ox;

V_p_f_f = M_f*(1 + 1/(B_f^(1/k_f) - 1))/rho_f;
V_p_ox_f = M_ox*(1 + 1/(B_ox^(1/k_ox) - 1))/rho_ox;

V_p_f_i = M_f/(rho_f*(B_f^(1/k_f) - 1));
V_p_ox_i = M_ox/(rho_ox*(B_ox^(1/k_ox) - 1));

m_f_i = m_p_i/(1 + OF_i);
m_ox_i = m_p_i*OF_i/(1 + OF_i);
dp_inj = alpha*p_c_i;
A_inj_f_tot = m_f_i/(C_d*sqrt(2*dp_inj*rho_f));
A_inj_ox_tot = m_ox_i/(C_d*sqrt(2*dp_inj*rho_ox));

L = h - L_c;
L_feed_f = 2*L/3;   % Assunto
L_feed_ox = L/3;    % Assunto
A_feed_f = pi*d_feed_f^2/4;
A_feed_ox = pi*d_feed_ox^2/4;
u_feed_f_i = m_f_i/(rho_f*A_feed_f);
u_feed_ox_i = m_ox_i/(rho_ox*A_feed_ox);

Re_f = rho_f*u_feed_f_i*d_feed_f/mu_f;
Re_ox = rho_ox*u_feed_ox_i*d_feed_ox/mu_ox;
f_f = moody(Re_f);
f_ox = moody(Re_ox);

K_f =  1 + f_f*L_feed_f/d_feed_f + (A_feed_f/(A_inj_f_tot*C_d))^2;
K_ox = 1 + f_ox*L_feed_ox/d_feed_ox + (A_feed_ox/(A_inj_ox_tot*C_d))^2;

dp_f = 0.5*rho_f*u_feed_f_i^2*K_f;
dp_ox = 0.5*rho_ox*u_feed_ox_i^2*K_ox;

p_f_i = p_c_i + dp_f;
p_ox_i = p_c_i + dp_ox;

% Iterazione
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
p_c = [p_c_i nan(1, length(tvet) - 1)];
c_t = [c_t_i nan(1, length(tvet) - 1)];
c_star = [c_star_i nan(1, length(tvet) - 1)];
T = [T_i*lambda nan(1, length(tvet) - 1)];
I_sp = [I_sp_i nan(1, length(tvet) - 1)];
T_f = [T_f_i nan(1, length(tvet) - 1)];
T_ox = [T_ox_i nan(1, length(tvet) - 1)];

% modifica Ale: inizializzo le variabili che mi servono 
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
E_rate=0.00001*10^-3 %erosion rate [mm/s]
m_p= [m_p_i nan(1, length(tvet) - 1)];
A_t_sal= [A_t nan(1, length(tvet) - 1)];
A_t_relativa = [A_t_relativa_i nan(1, length(tvet) - 1)];

V_f = 0;
V_ox = 0;
valido = 1;
valido_2 = 1;
cont = 1;
j = 1; % p_c_it inferiore di bisezione
toll = 0.1;

while valido
    % Modifiche Ale: considero ora la troat erosion e faccio cambiare At
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
    if V_f + V_ox + V_p_ox_i + V_p_f_i > V_tank_tot
        valido = 0;
        disp("Termine per volume occupato massimo raggiunto")
        continue
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

    while valido_2

        p_c_new = (p_c_up + p_c_dw)/2;

        u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
        u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
        
        m_f_new = rho_f*A_feed_f*u_feed_f_new;
        m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
        OF_new = m_ox_new/m_f_new;
        
        output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        c_star_cea = output.froz.cstar(1);
        c_t_new = output.froz.cf(end);
        T_c_new = output.froz.temperature(1);
        gamma_new = output.froz.gamma(1);
        c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
        T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star_new;
        I_sp_new = output.froz.isp(end);
        err = abs(c_star_cea - c_star_new);
        
        if abs(err) < toll
            valido_2 = 0;
        elseif  c_star_new < c_star_cea
            p_c_dw = p_c_new;
        else
            p_c_up = p_c_new;
        end
    end
    
    % Modifiche Ale: viscocity displaiment
    gamma_t_new=output.froz.gamma(2);
    rho_t_new=output.froz.density(2);
    mhu_t_new=output.froz.viscosity(2)*10^-4; %viscosità  
    mm_t=output.froz.mw(2); %massa molora propelente in gola 
    a_t_new=sqrt(gamma_t_new*output.froz.temperature(2)*(8.314/mm_t)); %velocità del suono alla throat  
    r_c_new=1.5*r_t(cont+1);   %curvatura alla gola  
    Re_t_new= rho_t_new*r_t(cont+1)*2*a_t_new/mhu_t_new; %reynolds in gola 
    Re_mod_new=sqrt((r_t(cont+1))/r_c_new)*Re_t_new; %Reynolds modificato 
    %discharge coefficient 
    Cd_new=1-(((gamma_t_new+1)/2)^(3/4))*(3.266-(2.128/(gamma_t_new+1)))*(Re_mod_new^(-1/2))+0.9428*((gamma_t_new-1)*(gamma_t_new+2)/((gamma_t_new+1)^(1/2)))*(Re_mod_new^-1); 
    % mi trovo m_p con m fuel + m ox
    m_re_new=Cd_new*(m_f_new + m_ox_new);
    A_t_relativa_new=(c_star_cea*m_re_new)/(output.froz.pressure(1)*10^5);

    valido_2 = 1;

    if p_c_new < p_c_min
        valido = 0;
        disp("Terminato per pressione in camera troppo bassa")
        continue
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
    c_star(cont + 1) = c_star_new;
    T(cont + 1) = T_new;
    I_sp(cont + 1) = I_sp_new;
    u_feed_f(cont + 1) = u_feed_f_new;
    u_feed_ox(cont + 1) = u_feed_ox_new;
    T_f(cont + 1) = T_f_new;
    T_ox(cont + 1) = T_ox_new;
    gamma(cont + 1) = gamma_new;

    %Modifica Ale: salvo le variabili 
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


    cont = cont + 1
end

% figure
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
title("Portate")
legend("m_f", "m_{ox}")

figure
plot(u_feed_f, 'r')
hold on
grid minor
plot(u_feed_ox, 'b')
title("Velocità feed")
legend("u_{feed,f}", "u_{feed,ox}")

figure
plot(I_sp)
grid minor
title("Impulso specifico")

figure
plot(T_f, 'r')
hold on
grid minor
plot(T_ox, 'b')
title("Temperature")
legend("Fuel", "Ossidante")

figure
plot(T_c)
title("Temperatura combustione")
grid minor

figure 
plot(Cd_Ale)
title("Variazione Cd")
grid minor

figure 
plot(m_p)
grid minor
hold on 
plot(m_re)
title("Mp e Mre")
legend("Mp","Mre")

figure 
plot(r_t)
grid minor
title("Variazione rt")

figure 
plot(A_t_sal)
grid minor
title("Varazione At")

figure 
plot(A_t_sal)
hold on
plot(A_t_relativa)
grid minor 
title("Confronto At reale e At relativa")
legend("At reale", "At relativa")


% T = T(~isnan(T));
% I_tot = sum(T)*dt

%(V_f + V_ox + V_N_f_i + V_He_ox_i)/V_tank_tot*100
