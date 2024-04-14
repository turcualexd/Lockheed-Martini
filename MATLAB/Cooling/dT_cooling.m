function [dT, q_flow] = dT_cooling(of, p_c, D_t, axial_distance, radius_vec,l_con, c_star, m_dot_fuel,c_RP1)

gamma_vec = [];
mach_vec = [];
Pr_CEA_vec = [];

sigma = @(T_rapp,M,gamma) 1./ ( (0.5*T_rapp.*(1 + M.^2 .* (gamma-1)/2) + 1/2).^(0.68) .* (1 + M.^2 .* (gamma - 1)/2 ).^0.12 );
n = length(radius_vec);
eps_vec = (radius_vec.^2 ./ (0.25*D_t^2));

x_cnv = linspace(0,l_con,10)';
y_cnv = D_t/2 + tan(pi/3).*(l_con - x_cnv);
axial_distance = axial_distance + l_con;
eps_c = (y_cnv/(D_t/2)).^2;
eps_ = flip(eps_c);
axial_distance = [x_cnv; axial_distance];
radius_vec = [y_cnv;radius_vec];

% convergent part
for j = 1:10
    eps_c(j) = (y_cnv(j)/(D_t/2))^2;
    output_c = cea(CEA('problem','rkt','eql','o/f',of,'fac','acat',eps_c,'case', 'DRY1', ...
        'p,bar',p_c/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423, 'wt%',100, ...
        'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
    gamma_vec  = [gamma_vec; output_c.eql.gamma(3:end)];
    Pr_CEA_vec = [Pr_CEA_vec; output_c.eql.prandtl.froz(3:end)];
    mach_vec   = [mach_vec; output_c.eql.mach(3:end)];

end

%divergent part
for i = 1:n
    eps = eps_vec(i);
    output = cea(CEA('problem','rkt','nfz',2,'o/f',of,'sup',eps,'case', 'DRY1', ...
        'p,bar',p_c/1e5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423, 'wt%',100, ...
        'oxid','O2(L)','O',2,'wt%',100,'output','massf','transport','trace',1e-10,'end'));
    if i ~= 1
        gamma_vec  = [gamma_vec; output.froz.gamma(3:end)];
        Pr_CEA_vec = [Pr_CEA_vec; output.froz.prandtl.froz(3:end)];
        mach_vec   = [mach_vec; output.froz.mach(3:end)];
    else
        gamma_vec  = [gamma_vec; output.froz.gamma(1:end-1)];
        Pr_CEA_vec = [Pr_CEA_vec; output.froz.prandtl.froz(1:end-1)];
        mach_vec   = [mach_vec; output.froz.mach(1:end-1)];
        T_c = output.froz.temperature(1);
        mu = output.froz.viscosity(1);
        c_p = output.froz.cp(1)*1e3; % J/Kg K
    end
end

r = Pr_CEA_vec.^1/3;
T_aw = T_c*( (1 + (mach_vec.^2).*r.*(gamma_vec - 1)./2) ./ (1 + (mach_vec.^2).*(gamma_vec - 1)./2));
T_wg_est = 1500; %estimated in order to see if the fuel is
% enough to absorb the power at reasonable termperature
T_wg_vec = (T_wg_est/T_aw(11)) * T_aw;
sigma_vec = sigma(T_wg_vec./T_c, mach_vec, gamma_vec);

% calculate heat_flux
for k = 1:length(eps_vec)
    A_conv = eps_c * (pi*0.25*D_t^2);
    A_div = eps_vec * (pi*0.25*D_t^2);
    A = [A_conv; A_div];
    h_g = h_g_BARTZ(mu, c_p, Pr_CEA_vec(1), D_t, p_c, c_star, 0.382*D_t/2,sigma_vec(2:end),A);
    q_dot_vec = h_g.*(T_aw(2:end) - T_wg_vec(2:end));
end
 q_flow = q_dot_vec;
%plot heat flux as a function of time

%calculate heat power along nozzle
Q = 0;
Q_L = 2*pi*q_dot_vec.*radius_vec;
for i = 1:length(axial_distance) - 1
    Q = Q + (Q_L(i+1) + Q_L(i))*(axial_distance(i+1) - axial_distance(i))*0.5;
end

dT = Q / (c_RP1*m_dot_fuel);
plot(axial_distance, q_dot_vec,'LineWidth',1.5);
xlabel('Axial Nozzle distance [m]'); ylabel('Heat flux [W/m^2]');
end

