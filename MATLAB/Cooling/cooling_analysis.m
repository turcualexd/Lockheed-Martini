clear; close all; clc;

n = 300;
D_t = 1.45e-2;

%interpolazione
[~, ~, ~, ~, ~, x3, y3] = rao_nozzle(D_t/2, 300, 100);
x3 = x3';
y3 = y3';
xx = linspace(x3(1),x3(end),n)';
yy = interp1(x3,y3,xx);
eps_c = 10;
eps_vec = (yy.^2 ./ (0.25*D_t^2));
OF = 2.33;
p_c = 50e5;

%f = @(x) (1./x).*( (1 + x.^2 .*(gamma-1)/2)./ (1 + (gamma-1)./2) ).^((gamma+1)./(2*(gamma-1))) - eps_vec;
%T_rapp_vec = linspace(0,1,5);
%sigma_vec = zeros(5,n);

gamma_vec = [];
mach_vec = [];
Pr_CEA_vec = [];
%%

for i = 1:n
    eps = eps_vec(i);
    output = cea(CEA('problem','rkt','nfz',2,'o/f',OF,'sup',eps,'case', 'DRY1', ...
        'p,bar',p_c/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423, 'wt%',100, ...
        'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
    if i ~= 1
        gamma_vec  = [gamma_vec; output.froz.gamma(3:end)];
        Pr_CEA_vec = [Pr_CEA_vec; output.froz.prandtl.froz(3:end)];
        mach_vec   = [mach_vec; output.froz.mach(3:end)];
        c_star = output.froz.cstar(end);
    else
        gamma_vec  = [gamma_vec; output.froz.gamma(1:end-1)];
        Pr_CEA_vec = [Pr_CEA_vec; output.froz.prandtl.froz(1:end-1)];
        mach_vec   = [mach_vec; output.froz.mach(1:end-1)];
        T_c = output.froz.temperature(1);
        mu = output.froz.viscosity(1);
        c_p = output.froz.cp(1)*1e3; % J/Kg K
    end
end

sigma = @(T_rapp,M,gamma) 1./ ( (0.5*T_rapp.*(1 + M.^2 .* (gamma-1)/2) + 1/2).^(0.68) .* (1 + M.^2 .* (gamma - 1)/2 ).^0.12 );

%

%Pr_for    = 4.*gamma_vec ./ (9.*gamma_vec - 5);
r = Pr_CEA_vec.^1/3;
T_aw = T_c*( (1 + (mach_vec.^2).*r.*(gamma_vec - 1)./2) ./ (1 + (mach_vec.^2).*(gamma_vec - 1)./2));
plot(eps_vec, T_aw(2:end))
T_wg_est = 1500; %estimated in order to see if the fuel is
% enough to absorb the power at reasonable termperature
c_star_design = c_star;

T_wg_vec = (T_wg_est/T_aw(1)) * T_aw;
sigma_vec = sigma(T_wg_vec./T_c, mach_vec, gamma_vec);

for k = 1:length(eps_vec)
    A = eps_vec * (pi*0.25*D_t^2);
    h_g = h_g_BARTZ(mu, c_p, Pr_CEA_vec(1), D_t, p_c, c_star_design, 0.382*D_t/2,sigma_vec(2:end),A);
    q_dot_vec = h_g.*(T_aw(2:end) - T_wg_vec(2:end));
end

%%

figure;
plot(eps_vec, q_dot_vec./1e3, "LineWidth",1.5);
xlabel("epsilon"); ylabel("heat flux [kW]");


TURCU = 0;
billy = 2*pi*q_dot_vec.*yy;
for i = 1:length(xx) - 1
    TURCU = TURCU + (billy(i+1) + billy(i))*(xx(i+1) - xx(i))*0.5;
end
m_dot_fuel = 0.085; %kg/s;
c  = 1880;
DT_fuel = TURCU / (c*m_dot_fuel);

%%
r_cc = (sqrt(eps_c)*D_t)/2;
l_con = (r_cc - D_t/2)/tan(pi/6);
dT = dT_cooling(OF, p_c, D_t, xx, yy, l_con, c_star, m_dot_fuel,c);
%notare  che questa condizione è all'equilibrio considerando però le
%condizioni iniziali. ricordarsi che il coefficiente h cala con la
%pressione in camera, quindi il flusso di calore diminuisce con il tempo
%%

 output_c = cea(CEA('problem','rkt','eql','o/f',2.24,'fac','acat',2 ,'case', 'DRY1', ...
        'p,bar',50,'reactants','fuel','RP-1(L)','C',1,'H',1.9423, 'wt%',100, ...
        'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));



