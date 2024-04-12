% clear, clc, close all
% 
% % Dati
% T_i = 1e3;          % N
% p_c_i = 50e5;       % Pa
% p_c_min = 20e5;     % Pa
% OF = 2.24;          % -
% T_f = 3571;         % K
% M_mol = 21.9e-3;    % kg/mol
% k = 1.24;           % -
% R_u = 8.3145;       % J/molK
% eps = 100;          % -
% g0 = 9.81;          % m/s^2
% l = 10000;
% p_vet = linspace(p_c_min, p_c_i, l);
% 
% % Calcoli
% R = R_u/M_mol;
% c1 = (2/(k + 1));
% c2 = (k + 1)/(k - 1);
% c3 = k/(k - 1);
% c4 = k*c3;
% c_star = sqrt(R*T_f/(k*c1^c2));
% m_p = zeros(1, l);
% I_s = m_p;
% for i = 1 : length(p_vet)
%     p_c = p_vet(i);
%     p_e = fzero(@(s) (1./c1).^(1./(k - 1)).*(s./p_c).^(1./k).*sqrt(c2.*(1 - (s./p_c).^(1./c3))) - 1./eps, 2700);
%     c_t = sqrt((2*c4*c1^c2)*(1 - (p_e/p_c)^(1/c3))) + eps*p_e/p_c;
%     A_t = T_i/(c_t*p_c);
%     A_e = eps*A_t;
%     u_e = sqrt(2*c3*R*T_f*(1 - (p_e/p_c)^(1/c3)));
%     m_p(i) = (T_i - p_e*A_e)/u_e;
%     I_s(i) = T_i/(m_p(i)*g0);
%     m_f = m_p(i)/(1 + OF);
%     m_ox = m_p(i)*OF/(1 + OF);
% end
% plot(p_vet, m_p)
% grid minor
% title("Portata")
% figure
% plot(p_vet, I_s)
% grid minor
% title("Impulso specifico")
% 
% %% Variamo O/F e p_c con CEAM
% clear, clc, close all
% 
% % Dati
% l = 100;
% OF_vet = linspace(1.5, 3, l+1);
% p_c_vet = linspace(50, 20, l-1);  % bar
% eps = 100;
% g0 = 9.81;
% T_i = 1000;     % N
% 
% % Calcoli
% gamma = zeros(length(OF_vet), length(p_c_vet));
% c_t = gamma;
% c_star = c_t;
% p_e = c_t;
% u_e = c_t;
% I_sp = c_t;
% m_p = c_t;
% T = c_t;
% %A_t = zeros(1, length(OF_vet));
% w = waitbar(0, 'CEAM goes brrrrrr...');
% cont = 1;
% for i = 1 : length(OF_vet)
%     for j = 1 : length(p_c_vet)
%         OF = OF_vet(i);
%         p_c = p_c_vet(j);
%         output = cea(CEA('problem','rkt','nfz',2,'o/f',OF,'sup',eps,'case','Porco Dio','p,bar',p_c,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
%         u_e(i, j) = output.froz.mach(end)*output.froz.sonvel(end);
%         p_e(i, j) = output.froz.pressure(end) * 1e5; % [Pa]
%         c_t(i, j) = output.froz.cf_vac(end);
%         c_star(i, j) = output.froz.cstar(end);
%         I_sp(i, j) = c_t(i, j)*c_star(i, j)/g0;
%         gamma(i, j) = output.froz.gamma(2);
%         if j == 1 && i == 1
%             m_p(i, j) = T_i/(c_t(i, j)*c_star(i, j));
%             A_t = m_p(i, j)/(output.froz.sonvel(2)*output.froz.density(2));
%             T(i, j) = T_i;
%         else
%             m_p(i, j) = A_t*p_c*1e5/c_star(i, j);
%             T(i, j) = m_p(i, j)*c_t(i, j)*c_star(i, j); 
%         end
%         waitbar(cont/l^2, w);
%         cont = cont + 1;
%     end
% end
% close(w)
% 
% %%
% clc, close all
% surf(p_c_vet, OF_vet, T);
% ylabel("O/F")
% xlabel("p_c")
% 
% %% 
% clear, clc, close all
% 
% % Dati iniziali noti
% T_i_n = 1000;       % N (Aggiungere perdite per ugello)
% p_c_i = 50e5;       % Pa
% p_c_min = 20e5;     % Pa
% rho_f = 807;        % kg/m^3
% rho_ox = 1140;      % kg/m^3
% L_star = 1.143;     % m
% h = 2;              % m
% d = 1;              % m
% g0 = 9.81;
% 
% % Dati assunti
% OF_i = 2.42;        % -
% eps = 300;          % -
% eps_c = 10;         % -
% A_max = 0.25*pi;    % m^2
% C_d = 0.7;          % -
% alpha = 0.2;        % -
% d_inj_f = 1e-3;     % m
% d_feed = 5e-3;      % m
% f = 0.02;           % -
% V_p_f_i = 0.2;     % m^3  
% V_p_ox_i = 0.2;    % m^3
% dt = 1;             % s
% lambda = 0.9974;
% k_He = 7/5;
% 
% % Dimensionamento a ritroso
% T_i = T_i_n/lambda;
% V_tot = pi*d^2*h/4;
% V_u = 0.8*V_tot;
% output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_i,'sup',eps,'case','Porco Dio','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
% c_star_i = output.froz.cstar(end);
% c_t_i = output.froz.cf(end);
% I_sp_i = c_t_i*c_star_i/g0;
% m_p_i = T_i/(c_t_i*c_star_i);
% A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
% A_e = eps*A_t;
% 
% if A_e > A_max
%     error("Area d'efflusso troppo grande")
% end
% 
% A_c = eps_c*A_t;
% d_c = 2*sqrt(A_c/pi);
% L_c = L_star/eps_c;
% V_int_c = pi*(d^2 - d_c^2)*L_c/4;
% 
% if V_int_c/V_tot > 0.2
%     h_tank = h - L_c;
% else
%     h_tank = h - L_c - 4*(0.2*V_tot - V_int_c)/(pi*d^2);
% end
% 
% V_tank_tot = pi*d^2*h_tank/4;
% m_f_i = m_p_i/(1 + OF_i);
% m_ox_i = m_p_i*OF_i/(1 + OF_i);
% dp_inj = alpha*p_c_i;
% A_f_tot = m_f_i/(C_d*sqrt(2*dp_inj*rho_f));
% A_ox_tot = m_ox_i/(C_d*sqrt(2*dp_inj*rho_ox));
% A_inj_f = pi*d_inj_f^2/4;
% N_inj_f = A_f_tot/A_inj_f;
% N_inj_f = floor(N_inj_f);
% N_inj_ox = 2*N_inj_f;
% A_inj_f = A_f_tot/N_inj_f;
% A_inj_ox = A_ox_tot/N_inj_ox;
% d_inj_f = 2*sqrt(A_inj_f/pi);
% d_inj_ox = 2*sqrt(A_inj_ox/pi);
% v_inj_f = m_f_i/(rho_f*A_f_tot);
% v_inj_ox = m_ox_i/(rho_ox*A_ox_tot);
% 
% L = h - L_c;
% L_feed_f = 2*L/3;   % Assunto
% L_feed_ox = L/3;    % Assunto
% A_feed = pi*d_feed^2/4;
% u_feed_f = m_f_i/(rho_f*A_feed);
% u_feed_ox = m_ox_i/(rho_ox*A_feed);
% dp_feed_f = 0.5*rho_f*f*u_feed_f^2*L_feed_f/d_feed;
% dp_feed_ox = 0.5*rho_ox*f*u_feed_ox^2*L_feed_ox/d_feed;
% dp_tank_f = 0.5*rho_f*u_feed_f^2;
% dp_tank_ox = 0.5*rho_ox*u_feed_ox^2;
% dp_f = dp_inj + dp_feed_f + dp_tank_f;
% dp_ox = dp_inj + dp_feed_ox + dp_tank_ox;
% p_p_f_i = p_c_i + dp_f;
% p_p_ox_i = p_c_i + dp_ox;
% 
% % Espansione
% p_p_f = p_p_f_i;
% p_p_ox = p_p_ox_i;
% V_p_f = V_p_f_i;
% V_p_ox = V_p_ox_i;
% m_f = m_f_i;
% m_ox = m_ox_i;
% t = 0;
% V_f = 0;
% V_ox = 0;
% p_c = p_c_i;
% OF = OF_i;
% c_t = c_t_i;
% c_star = c_star_i;
% T = T_i_n;
% I_sp = I_sp_i;
% %alpha_ox = alpha;
% valido = 1;
% 
% opt = optimoptions("fsolve", "FunctionTolerance", 1e-15);
% 
% while valido
%     V_f_dt = m_f(end)*dt/rho_f;
%     V_p_f_new = V_p_f(end) + V_f_dt;
%     V_f = V_f + V_f_dt;
%     p_p_f_new = p_p_f(end)/(1 + V_f_dt/V_p_f(end))^k_He;
%     V_ox_dt = m_ox(end)*dt/rho_ox;
%     V_p_ox_new = V_p_ox(end) + V_ox_dt; 
%     V_ox = V_ox + V_ox_dt;
%     V = V_f + V_ox + V_p_ox_i + V_p_f_i;
% 
%     if V > V_tank_tot
%         valido = 0;
%         disp("Termine per volume occupato massimo raggiunto")
%         continue
%     end
% 
%     p_p_ox_new = p_p_ox(end)/(1 + V_ox_dt/V_p_ox(end))^k_He;
%     fun = @(x) [alpha/(1 + alpha)*(p_p_f_new - 0.5*rho_f*x(1).^2*(1 + f*L_feed_f/d_feed)) - 0.5*rho_f*(A_feed*x(1)./(x(3)*A_f_tot));
%                 alpha/(1 + alpha)*(p_p_ox_new - 0.5*rho_ox*x(2).^2*(1 + f*L_feed_f/d_feed)) - 0.5*rho_ox*(A_feed*x(2)./(x(3)*A_ox_tot));
%                 (p_p_f_new - p_p_ox_new) + 0.5*(rho_ox*x(2).^2*(1 + f*L_feed_ox/d_feed) - rho_f*x(1).^2*(1 + f*L_feed_f/d_feed)) + 0.5*rho_ox*(A_feed*x(2)./(x(3)*A_ox_tot)).^2 - 0.5*rho_f*(A_feed*x(1)./(x(3)*A_f_tot)).^2];
%     x = fsolve(fun, [u_feed_f(end); u_feed_f(end); C_d(end)+0.5], opt);
%     u_feed_f_new = x(1);
%     u_feed_ox_new = x(2);
%     C_d_new = x(3);
%     % u_feed_f = sqrt((alpha*p_p_f_new/(1 + alpha))/(0.5*rho_f + 0.5*rho_f*f*L_feed_f/d_feed + 0.5*rho_f*(A_feed/(C_d*A_f_tot))^2));
%     dp_f = 0.5*rho_f*u_feed_f_new^2*(1 + f*L_feed_f/d_feed) + 0.5*rho_f*(u_feed_f_new*A_feed/(C_d_new*A_f_tot))^2;
%     p_c_new = p_p_f_new - dp_f;
% 
%     if p_c_new < p_c_min
%         valido = 0;
%         disp("Terminato per pressione in camera troppo bassa")
%         continue
%     end
% 
%     % u_feed_ox = sqrt((p_p_ox_new - p_c_new)/(0.5*rho_ox*(1 + f*L_feed_ox/d_feed + (A_feed/(C_d*A_ox_tot))^2)));
%     % dp_tfd_ox = 0.5*rho_ox*u_feed_ox^2*(1 + f*L_feed_ox/d_feed);
%     % beta = (p_p_ox_new - dp_tfd_ox - p_c_new)/(p_p_ox_new - dp_tfd_ox);
%     % alpha_ox_new = beta/(1 - beta);
%     m_f_new = rho_f*A_feed*u_feed_f_new;
%     m_ox_new = rho_ox*A_feed*u_feed_ox_new;
%     OF_new = m_ox_new/m_f_new;
% 
%     if mod(t, 20) == 0
%         output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
%         c_t_new = output.froz.cf(end);
%         c_star_new = output.froz.cstar(end);
%         T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star_new;
%         I_sp_new = c_t_new*c_star_new/g0;
% 
%         c_t = [c_t c_t_new];
%         c_star = [c_star c_star_new];
%         T = [T T_new];
%         I_sp = [I_sp I_sp_new];
%     end
% 
%     V_p_f = [V_p_f V_p_f_new];
%     V_p_ox = [V_p_ox V_p_ox_new];
%     p_p_f = [p_p_f p_p_f_new];
%     p_p_ox = [p_p_ox p_p_ox_new];
%     p_c = [p_c p_c_new];
%     %alpha_ox = [alpha_ox alpha_ox_new];
%     u_feed_f = [u_feed_f u_feed_f_new];
%     u_feed_ox = [u_feed_ox u_feed_ox_new];
%     C_d = [C_d C_d_new];
%     m_f = [m_f m_f_new];
%     m_ox = [m_ox m_ox_new];
%     OF = [OF OF_new];
%     t = t + dt;
% end
% 
% tvet = 0 : dt : t;
% plot(t, C_d, 'r')
% grid minor

%% Dimensionamento iterativo
% clear, clc, close all
% 
% % Dati iniziali noti
% T_i_n = 1000;       % N (Aggiungere perdite per ugello)
% p_c_i = 50e5;       % Pa
% p_c_min = 20e5;     % Pa
% rho_f = 807;        % kg/m^3
% rho_ox = 1140;      % kg/m^3
% L_star = 1.143;     % m
% h = 2;              % m
% d = 1;              % m
% g0 = 9.81;          % m/s^2
% R = 8.314;          % J/molK
% M_m_He = 4e-3;      % kg/mol
% M_m_N = 28e-3;      % kg/mol
% t_max = 10000;      % s
% mu_f = 0.75e-3;     % Pas
% mu_ox = 0.196e-3;   % Pas
% 
% % Dati assunti
% OF_i = 2.42;           % -
% eps = 300;          % -
% eps_c = 10;         % -
% A_max = 0.25*pi;    % m^2
% C_d = 0.7;          % -
% alpha = 0.2;        % -
% d_inj_f = 1e-3;     % m
% d_feed_f = 5e-3;    % m
% d_feed_ox = 7e-3;   % m
% dt = 1;             % s
% lambda = 0.9974;    % -
% k_He = 5/3;         % -
% k_N = 7/5;          % -
% T_N_f_i = 300;      % K
% T_He_ox_i = 90;     % K
% V_N_f_i = 0.2;      % m^3
% V_He_ox_i = 0.2;    % m^3
% 
% % Dimensionamento a ritroso
% T_i = T_i_n/lambda;
% V_tot = pi*d^2*h/4;
% output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_i,'sup',eps,'case','Porco Dio','p,bar',p_c_i/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
% c_star_i = output.froz.cstar(end);
% c_t_i = output.froz.cf(end);
% T_c_i = output.froz.temperature(1);
% gamma_i = output.froz.gamma(1);
% I_sp_i = output.froz.isp(end);
% m_p_i = T_i/(c_t_i*c_star_i);
% A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
% A_e = eps*A_t;
% 
% if A_e > A_max
%     error("Area d'efflusso troppo grande")
% end
% 
% A_c = eps_c*A_t;
% d_c = 2*sqrt(A_c/pi);
% L_c = L_star/eps_c;
% V_int_c = pi*(d^2 - d_c^2)*L_c/4;
% 
% if V_int_c/V_tot > 0.2
%     h_tank = h - L_c;
% else
%     h_tank = h - L_c - 4*(0.2*V_tot - V_int_c)/(pi*d^2);
% end
% 
% V_tank_tot = pi*d^2*h_tank/4;
% m_f_i = m_p_i/(1 + OF_i);
% m_ox_i = m_p_i*OF_i/(1 + OF_i);
% dp_inj = alpha*p_c_i;
% A_inj_f_tot = m_f_i/(C_d*sqrt(2*dp_inj*rho_f));
% A_inj_ox_tot = m_ox_i/(C_d*sqrt(2*dp_inj*rho_ox));
% 
% L = h - L_c;
% L_feed_f = 2*L/3;   % Assunto
% L_feed_ox = L/3;    % Assunto
% A_feed_f = pi*d_feed_f^2/4;
% A_feed_ox = pi*d_feed_ox^2/4;
% u_feed_f_i = m_f_i/(rho_f*A_feed_f);
% u_feed_ox_i = m_ox_i/(rho_ox*A_feed_ox);
% 
% Re_f = rho_f*u_feed_f_i*d_feed_f/mu_f;
% Re_ox = rho_ox*u_feed_ox_i*d_feed_ox/mu_ox;
% f_f = moody(Re_f);
% f_ox = moody(Re_ox);
% 
% K_f =  1 + f_f*L_feed_f/d_feed_f + (A_feed_f/(A_inj_f_tot*C_d))^2;
% K_ox = 1 + f_ox*L_feed_ox/d_feed_ox + (A_feed_ox/(A_inj_ox_tot*C_d))^2;
% 
% dp_f = 0.5*rho_f*u_feed_f_i^2*K_f;
% dp_ox = 0.5*rho_ox*u_feed_ox_i^2*K_ox;
% 
% p_f_i = p_c_i + dp_f;
% p_ox_i = p_c_i + dp_ox;
% 
% % Iterazione
% tvet = 0 : dt : t_max;
% m_f = [m_f_i nan(1, length(tvet) - 1)];
% m_ox = [m_ox_i nan(1, length(tvet) - 1)];
% u_feed_f = [u_feed_f_i nan(1, length(tvet) - 1)];
% u_feed_ox = [u_feed_ox_i nan(1, length(tvet) - 1)];
% V_N_f = [V_N_f_i nan(1, length(tvet) - 1)];
% V_He_ox = [V_He_ox_i nan(1, length(tvet) - 1)];
% p_f = [p_f_i nan(1, length(tvet) - 1)];
% p_ox = [p_ox_i nan(1, length(tvet) - 1)];
% p_c = [p_c_i nan(1, length(tvet) - 1)];
% OF = [OF_i nan(1, length(tvet) - 1)];
% T_c = [T_c_i nan(1, length(tvet) - 1)];
% gamma = [gamma_i nan(1, length(tvet) - 1)];
% p_c = [p_c_i nan(1, length(tvet) - 1)];
% c_t = [c_t_i nan(1, length(tvet) - 1)];
% c_star = [c_star_i nan(1, length(tvet) - 1)];
% T = [T_i*lambda nan(1, length(tvet) - 1)];
% I_sp = [I_sp_i nan(1, length(tvet) - 1)];
% T_N_f = [T_N_f_i nan(1, length(tvet) - 1)];
% T_He_ox = [T_He_ox_i nan(1, length(tvet) - 1)];
% 
% V_f = 0;
% V_ox = 0;
% valido = 1;
% valido_2 = 1;
% cont = 1;
% j = 1; % p_c_it inferiore di bisezione
% toll = 0.1;
% 
% while valido
%     dV_f = m_f(cont)*dt/rho_f;
%     dV_ox = m_ox(cont)*dt/rho_ox;
% 
%     V_N_f_new = V_N_f(cont) + dV_f;
%     V_He_ox_new = V_He_ox(cont) + dV_ox;
% 
%     T_N_f_new = T_N_f(cont)*(V_N_f(cont)/V_N_f_new)^(k_N - 1);
%     T_He_ox_new = T_He_ox(cont)*(V_He_ox(cont)/V_He_ox_new)^(k_He - 1);
% 
%     p_f_new = p_f(cont)*(V_N_f(cont)/V_N_f_new)^k_N;
%     p_ox_new = p_ox(cont)*(V_He_ox(cont)/V_He_ox_new)^k_He;
% 
%     V_f = V_f + dV_f;
%     V_ox = V_ox + dV_ox;
%     if V_f + V_ox + V_He_ox_i + V_N_f_i > V_tank_tot
%         valido = 0;
%         disp("Termine per volume occupato massimo raggiunto")
%         continue
%     end
% 
%     p_c_new = p_c(cont);
%     u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
%     u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
% 
%     m_f_new = rho_f*A_feed_f*u_feed_f_new;
%     m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
%     OF_new = m_ox_new/m_f_new;
% 
%     output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
%     c_star_cea = output.froz.cstar(1);
%     c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
% 
%     while c_star_new > c_star_cea 
% 
%         j = j - 0.01;
%         p_c_new = j*p_c(cont);
% 
%         u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
%         u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
% 
%         m_f_new = rho_f*A_feed_f*u_feed_f_new;
%         m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
%         OF_new = m_ox_new/m_f_new;
% 
%         output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
%         c_star_cea = output.froz.cstar(1);
%         c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
% 
%     end
% 
%     p_c_up = (j + 0.01)*p_c(cont);
%     p_c_dw = j*p_c(cont);
%     j = 1;
% 
%     while valido_2
% 
%         p_c_new = (p_c_up + p_c_dw)/2;
% 
%         u_feed_f_new = sqrt(2*(p_f_new - p_c_new)/(rho_f*K_f));
%         u_feed_ox_new = sqrt(2*(p_ox_new - p_c_new)/(rho_ox*K_ox));
% 
%         m_f_new = rho_f*A_feed_f*u_feed_f_new;
%         m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
%         OF_new = m_ox_new/m_f_new;
% 
%         output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
%         c_star_cea = output.froz.cstar(1);
%         c_t_new = output.froz.cf(end);
%         T_c_new = output.froz.temperature(1);
%         gamma_new = output.froz.gamma(1);
%         c_star_new = A_t*p_c_new/(m_f_new + m_ox_new);
%         T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star_new;
%         I_sp_new = output.froz.isp(end);
%         err = abs(c_star_cea - c_star_new);
% 
%         if abs(err) < toll
%             valido_2 = 0;
%         elseif  c_star_new < c_star_cea
%             p_c_dw = p_c_new;
%         else
%             p_c_up = p_c_new;
%         end
%     end
% 
%     valido_2 = 1;
% 
%     if p_c_new < p_c_min
%         valido = 0;
%         disp("Terminato per pressione in camera troppo bassa")
%         continue
%     end
% 
%     V_N_f(cont + 1) = V_N_f_new;
%     V_He_ox(cont + 1) = V_He_ox_new;
%     p_f(cont + 1) = p_f_new;
%     p_ox(cont + 1) = p_ox_new;
%     m_f(cont + 1) = m_f_new;
%     m_ox(cont + 1) = m_ox_new;
%     OF(cont + 1) = OF_new;
%     T_c(cont + 1) = T_c_new;
%     p_c(cont + 1) = p_c_new;
%     c_t(cont + 1) = c_t_new;
%     c_star(cont + 1) = c_star_new;
%     T(cont + 1) = T_new;
%     I_sp(cont + 1) = I_sp_new;
%     u_feed_f(cont + 1) = u_feed_f_new;
%     u_feed_ox(cont + 1) = u_feed_ox_new;
%     T_N_f(cont + 1) = T_N_f_new;
%     T_He_ox(cont + 1) = T_He_ox_new;
%     gamma(cont + 1) = gamma_new;
% 
%     cont = cont + 1
% end
% 
% % figure
% plot(OF)
% grid minor
% title("OF")
% 
% figure
% plot(p_f, 'r')
% grid minor
% hold on
% plot(p_ox, 'b')
% plot(p_c, 'm')
% title("Pressioni")
% legend("p_f", "p_{ox}", "p_c")
% 
% figure
% plot(m_f, 'r')
% hold on
% grid minor
% plot(m_ox, 'b')
% title("Portate")
% legend("m_f", "m_{ox}")
% 
% figure
% plot(u_feed_f, 'r')
% hold on
% grid minor
% plot(u_feed_ox, 'b')
% title("Velocità feed")
% legend("u_{feed,f}", "u_{feed,ox}")
% 
% figure
% plot(I_sp)
% grid minor
% title("Impulso specifico")
% 
% figure
% plot(T_N_f, 'r')
% hold on
% grid minor
% plot(T_He_ox, 'b')
% title("Temperature")
% legend("Fuel", "Ossidante")
% 
% figure
% plot(gamma)
% title("gamma comb")
% grid minor
% 
% % T = T(~isnan(T));
% % I_tot = sum(T)*dt
% 
% (V_f + V_ox + V_N_f_i + V_He_ox_i)/V_tank_tot*100

%% C_star costante
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
c_star = output.froz.cstar(end);
c_t_i = output.froz.cf(end);
T_c_i = output.froz.temperature(1);
gamma_i = output.froz.gamma(1);
I_sp_i = output.froz.isp(end);
m_p_i = T_i/(c_t_i*c_star);
A_t = m_p_i/(output.froz.sonvel(2)*output.froz.density(2));
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
T = [T_i*lambda nan(1, length(tvet) - 1)];
I_sp = [I_sp_i nan(1, length(tvet) - 1)];
T_f = [T_f_i nan(1, length(tvet) - 1)];
T_ox = [T_ox_i nan(1, length(tvet) - 1)];

V_f = 0;
V_ox = 0;
valido = 1;
cont = 1;
j = 1; % p_c_it inferiore di bisezione

while valido
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

    fun = @(x) [x(3) - p_f_new + 0.5*rho_f*K_f*x(1)^2;
                x(3) - p_ox_new + 0.5*rho_ox*K_ox*x(2)^2;
                c_star*(rho_f*A_feed_f*x(1) + rho_ox*A_feed_ox*x(2)) - A_t*x(3)];
    x = fsolve(fun, [u_feed_f(cont), u_feed_ox(cont), p_c(cont)], optimoptions("fsolve", "Display", "none"));
    u_feed_f_new = x(1);
    u_feed_ox_new = x(2);
    p_c_new = x(3);

    if p_c_new < p_c_min
        valido = 0;
        disp("Terminato per pressione in camera troppo bassa")
        continue
    end
    
    m_f_new = rho_f*A_feed_f*u_feed_f_new;
    m_ox_new = rho_ox*A_feed_ox*u_feed_ox_new;
    OF_new = m_ox_new/m_f_new;
    
    output = cea(CEA('problem','rkt','nfz',2,'o/f',OF_new,'sup',eps,'case','Porco Dio','p,bar',p_c_new/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
    c_t_new = output.froz.cf(end);
    T_c_new = output.froz.temperature(1);
    I_sp_new = output.froz.isp(end);
    gamma_new = output.froz.gamma(1);
    T_new = lambda*(m_f_new + m_ox_new)*c_t_new*c_star;

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
    T(cont + 1) = T_new;
    I_sp(cont + 1) = I_sp_new;
    u_feed_f(cont + 1) = u_feed_f_new;
    u_feed_ox(cont + 1) = u_feed_ox_new;
    T_f(cont + 1) = T_f_new;
    T_ox(cont + 1) = T_ox_new;
    gamma(cont + 1) = gamma_new;

    cont = cont + 1
end

figure
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
plot(gamma)
title("gamma comb")
grid minor

T = T(~isnan(T));
I_tot = sum(T)*dt;
(V_f + V_ox + V_N_f_i + V_He_ox_i)/V_tank_tot*100


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


% T = T(~isnan(T));
% I_tot = sum(T)*dt

%(V_f + V_ox + V_N_f_i + V_He_ox_i)/V_tank_tot*100
