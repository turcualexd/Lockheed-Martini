clear, clc, close all


% questo ha i VAC corretti
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
dt = 120;             % s 
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
V_loss = 0.25*pi*((L_c + L_con)*d^2 - L_c*d_c^2 - L_con(d_c^2 + d_c*d_t + d_t^2)/3);

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

%qua incertezze
coeff = 96.7/1903.3;

sigma_f = coeff * d_inj_f;
sigma_ox = coeff * d_inj_ox;

n_simulations = 15; %da aumentare dopo

tvet = 0 : dt : t_max;
m_f_mat = [nan(n_simulations, length(tvet))];
m_ox_mat = [nan(n_simulations, length(tvet))];
u_feed_f_mat = [nan(n_simulations, length(tvet))];
u_feed_ox_mat = [nan(n_simulations, length(tvet))];
V_p_f_mat = [nan(n_simulations, length(tvet))];
V_p_ox_mat = [nan(n_simulations, length(tvet))];
p_f_mat = [nan(n_simulations, length(tvet))];
p_ox_mat = [nan(n_simulations, length(tvet))];
p_c_mat = [nan(n_simulations, length(tvet))];
OF_mat = [nan(n_simulations, length(tvet))];
T_c_mat = [nan(n_simulations, length(tvet))];
gamma_mat = [nan(n_simulations, length(tvet))];
c_t_mat = [nan(n_simulations, length(tvet))];
c_star_mat = [nan(n_simulations, length(tvet))];
T_mat = [nan(n_simulations, length(tvet))];
I_sp_mat = [nan(n_simulations, length(tvet))];
T_f_mat = [nan(n_simulations, length(tvet))];
T_ox_mat = [nan(n_simulations, length(tvet))];

w = waitbar(0, 'CEAM goes brrrrrr... (0%)');
contatore = 1;

for q = 1 : n_simulations
    d_inj_f_vec = nan(N_inj_f,1);
    d_inj_ox_vec = nan(N_inj_ox,1);
    A_inj_f_vec = nan(N_inj_f,1);
    A_inj_ox_vec = nan(N_inj_ox,1);

    % ANDRA CICLATO PER DIVERSE ITERAZIONI PER AVERE PIU POSSIBILITA DI VARIAZIONI
    % E ANDRA CICLATO PER DIVERSI CD


    for i = 1 : N_inj_f
        d_inj_f_vec(i) = normrnd(d_inj_f, sigma_f);
        A_inj_f_vec(i) = pi*(d_inj_f_vec(i))^2/4;
    end
    for i = 1 : N_inj_ox
        d_inj_ox_vec(i) = normrnd(d_inj_ox, sigma_ox);
        A_inj_ox_vec(i) = pi*(d_inj_ox_vec(i))^2/4;
    end


    % sovrascrivo area totale che diventa diversa
    A_inj_f_tot = sum(A_inj_f_vec);
    A_inj_ox_tot = sum(A_inj_ox_vec);

    % sovrascrivo con nuove portate di f e ox

    m_f_i = A_inj_f_tot * (C_d*sqrt(2*dp_inj*rho_f));
    m_ox_i = A_inj_ox_tot * (C_d*sqrt(2*dp_inj*rho_ox));

    OF_i = m_ox_i/m_f_i;

    m_p_i = m_f_i + m_ox_i;



    %%

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
        T(cont + 1) = T_new;
        I_sp(cont + 1) = I_sp_new;
        u_feed_f(cont + 1) = u_feed_f_new;
        u_feed_ox(cont + 1) = u_feed_ox_new;
        T_f(cont + 1) = T_f_new;
        T_ox(cont + 1) = T_ox_new;
        gamma(cont + 1) = gamma_new;

        cont = cont + 1
    end
    m_f_mat(q, :) = m_f;
    m_ox_mat(q, :) = m_ox;
    u_feed_f_mat(q, :) = u_feed_f;
    u_feed_ox_mat(q, :) = u_feed_ox;
    V_p_f_mat(q, :) = V_p_f;
    V_p_ox_mat(q, :) = V_p_ox;
    p_f_mat(q, :) = p_f;
    p_ox_mat(q, :) = p_ox;
    p_c_mat(q, :) = p_c;
    OF_mat(q, :) = OF;
    T_c_mat(q, :) = T_c;
    gamma_mat(q, :) = gamma;
    c_t_mat(q, :) = c_t;
    c_star_mat(q, :) = c_star;
    T_mat(q, :) = T;
    I_sp_mat(q, :) = I_sp;
    T_f_mat(q, :) = T_f;
    T_ox_mat(q, :) = T_ox;

    waitbar(contatore/n_simulations, w, sprintf('CEAM goes brrrrrr... (%.2f%%)', 100*contatore/n_simulations));
    contatore = contatore + 1;
end
close(w)

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

%%
figure
hold on
grid minor
plot(tvet,OF_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, OF_mat(q, :), 'Color','k');
end
title("O/F Ratio")
xlabel("t [s]")
ylabel("O/F [-]")
legend('avg', 'sim')
%%
figure
hold on
grid minor
plot(tvet,p_ox_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, p_ox_mat(q, :), 'Color','k');
end
title("Oxidizer pressure")
xlabel("t [s]")
ylabel("p_{ox} [Pa]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,p_f_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, p_f_mat(q, :), 'Color','k');
end
title("Fuel pressure")
xlabel("t [s]")
ylabel("p_{f} [Pa]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,p_c_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, p_c_mat(q, :), 'Color','k');
end
title("Chamber pressure")
xlabel("t [s]")
ylabel("p_{c} [Pa]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,m_ox_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, m_ox_mat(q, :), 'Color','k');
end
title("Oxidizer mass flow rate")
xlabel("t [s]")
ylabel("m_{ox} [kg/s]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,m_f_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, m_f_mat(q, :), 'Color','k');
end
title("Fuel mass flow rate")
xlabel("t [s]")
ylabel("m_{f} [kg/s]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,u_feed_ox_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, u_feed_ox_mat(q, :), 'Color','k');
end
title("Oxidizer feed velocity")
xlabel("t [s]")
ylabel("u_{feed,ox} [m/s]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,u_feed_f_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, u_feed_f_mat(q, :), 'Color','k');
end
title("Fuel feed velocity")
xlabel("t [s]")
ylabel("u_{feed,f} [m/s]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,I_sp_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, I_sp_mat(q, :), 'Color','k');
end
title("Specific Impulse")
xlabel("t [s]")
ylabel("I_{sp} [s]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,T_ox_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, T_ox_mat(q, :), 'Color','k');
end
title("Oxidizer temperature")
xlabel("t [s]")
ylabel("T_{oxidizer} [K]")
legend('avg', 'sim')

figure
hold on
grid minor
plot(tvet,T_f_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, T_f_mat(q, :), 'Color','k');
end
title("Fuel temperature")
xlabel("t [s]")
ylabel("T_{fuel} [K]")
legend('avg', 'sim')


figure
hold on
grid minor
plot(tvet,T_c_avg, 'LineWidth',2, 'Color','r')
for q = 1:n_simulations
    plot(tvet, T_c_mat(q, :), 'Color','k');
end
title("Combustion Chamber temperature")
xlabel("t [s]")
ylabel("T_{cc} [K]")
legend('avg', 'sim')



dp_c_end = p_c_new - p_c_min;
I_tot = sum(T)*dt;
V_occ = (V_f - dV_f + V_ox - dV_ox + V_p_f_i + V_p_ox_i)/V_tank_tot*100;

R_p_f = R/M_m_p_f;
R_p_ox = R/M_m_p_ox;

M_p_f = p_f_i*V_p_f_i/(R_p_f*T_f_i);
M_p_ox = p_ox_i*V_p_ox_i/(R_p_ox*T_ox_i);