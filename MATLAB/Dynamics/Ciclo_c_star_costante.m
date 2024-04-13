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
t_max = 4000;       % s
mu_f = 0.75e-3;     % Pas
mu_ox = 0.196e-3;   % Pas

% Dati assunti
OF_vet = 2.4 : 0.01 : 2.5;
eps = 300;          % -
eps_c = 10;         % -
C_d = 0.82;         % -
alpha = 0.2;        % -
d_feed_f = 5e-3;    % m
d_feed_ox = 7e-3;   % m
dt = 1;             % s
lambda = 1;         % -
k_ox = 7/5;         % -
k_f = 7/5;          % -
B_vet = 2.79;       % -
alpha_con = 30;     % deg

% Dimensionamento a ritroso
T_i = T_i_n/lambda;

V_tot = pi*d^2*h/4;

w = waitbar(0, 'CEAM goes brrrrrr... (0%)');
contatore = 1;
I_tot = nan(length(OF_vet), length(B_vet));

for i = 1 : length(OF_vet)
    for j = 1 : length(B_vet)
        
        OF_i = OF_vet(i);
        OF_m = OF_i;
        
        B_f = B_vet(j);
        B_ox = B_f;
        
        output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_i,'sup',eps,'case','Porco Dio','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
        
        c_star = output.froz.cstar(end);
        c_t_i = output.froz.cf(end);
        m_p_i = T_i/(c_t_i*c_star);
        
        A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
        A_c = eps_c*A_t;
        
        d_t = 2*sqrt(A_t/pi);
        d_c = 2*sqrt(A_c/pi);
        
        L_con = (d_c - d_t)/(2*tand(alpha_con));
        
        L_c = L_star/eps_c;
        V_loss = 0.25*pi*((L_c + L_con)*d^2 - L_c*d_c^2 - L_con*(d_c^2 + d_c*d_t + d_t^2)/3);
        
        if V_loss/V_tot > 0.2
            h_t = h - L_c - L_con;
        else
            h_t = h - L_c - L_con - 4*(0.2*V_tot - V_loss)/(pi*d^2);
        end
        
        V_tank_tot = pi*d^2*h_t/4;

        % Sistema Alex
        M_f = V_tank_tot/(OF_m*(1 + 1/(B_ox^(1/k_ox) - 1))/rho_ox + (1 + 1/(B_f^(1/k_f) - 1))/rho_f);
        M_ox = OF_m*M_f;
        
        V_t_f = M_f*(1 + 1/(B_f^(1/k_f) - 1))/rho_f;
        h_t_f = 4*V_t_f/(pi*d^2);
        
        L_feed_f = h - L_c - h_t_f;
        L_feed_ox = h - L_c - h_t;    

        V_p_f_i = M_f/(rho_f*(B_f^(1/k_f) - 1));
        V_p_ox_i = M_ox/(rho_ox*(B_ox^(1/k_ox) - 1));
        
        m_f_i = m_p_i/(1 + OF_i);
        m_ox_i = m_p_i*OF_i/(1 + OF_i);
        
        dp_inj = alpha*p_c_i;
        
        A_inj_f_tot = m_f_i/(C_d*sqrt(2*dp_inj*rho_f));
        A_inj_ox_tot = m_ox_i/(C_d*sqrt(2*dp_inj*rho_ox));
        
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
        T = [T_i*lambda nan(1, length(tvet) - 1)];
        
        V_f = 0;
        V_ox = 0;
        valido = 1;
        cont = 1;
        
        while true
            dV_f = m_f(cont)*dt/rho_f;
            dV_ox = m_ox(cont)*dt/rho_ox;
        
            V_p_f_new = V_p_f(cont) + dV_f;
            V_p_ox_new = V_p_ox(cont) + dV_ox;
        
            p_f_new = p_f(cont)*(V_p_f(cont)/V_p_f_new)^k_f;
            p_ox_new = p_ox(cont)*(V_p_ox(cont)/V_p_ox_new)^k_ox;
        
            V_f = V_f + dV_f;
            V_ox = V_ox + dV_ox;
            
            if V_f > M_f/rho_f || V_ox > M_ox/rho_ox
                break
            end
        
            fun = @(x) c_star - A_t*x/(A_feed_f*sqrt(2*rho_f*(p_f_new - x)/K_f) + A_feed_ox*sqrt(2*rho_ox*(p_ox_new - x)/K_ox));
            p_c_new = fsolve(fun, p_c(cont), optimoptions("fsolve", "Display", "none"));
            
            if p_c_new < p_c_min
                break
            end

            u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
            u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
            
            m_f_new = rho_f*A_feed_f*u_feed_f_new;
            m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
            OF_new = m_ox_new/m_f_new;
            
            output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
            c_t_new = output.froz.cf(end);
            T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star;
        
            V_p_f(cont + 1) = V_p_f_new;
            V_p_ox(cont + 1) = V_p_ox_new;
            p_f(cont + 1) = p_f_new;
            p_ox(cont + 1) = p_ox_new;
            m_f(cont + 1) = m_f_new;
            m_ox(cont + 1) = m_ox_new;
            p_c(cont + 1) = p_c_new;
            T(cont + 1) = T_new;
            u_feed_f(cont + 1) = u_feed_f_new;
            u_feed_ox(cont + 1) = u_feed_ox_new;
        
            cont = cont + 1;
        end
        T = T(~isnan(T));        
        I_tot(i, j) = sum(T)*dt;

        waitbar(contatore/(length(OF_vet)*length(B_vet)), w, sprintf('CEAM goes brrrrrr... (%.2f%%)', 100*contatore/(length(OF_vet)*length(B_vet))));
        contatore = contatore + 1;
    end
end
close(w)

[m, i] = max(I_tot);
[I_tot_best, j] = max(m);
i = i(j);

OF_best = OF_vet(i)
B_best = B_vet(j)
I_tot_best