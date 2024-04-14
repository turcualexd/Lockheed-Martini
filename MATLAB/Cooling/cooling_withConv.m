clear; close all; clc;

n = 300;
D_t = 1.45e-2;

%interpolazione
[~, ~, ~, ~, ~, x3, y3] = rao_nozzle(D_t/2, 300, 100);
x3 = x3';
y3 = y3';
xx = linspace(x3(1),x3(end),n)';
yy = interp1(x3,y3,xx);
eps_vec = (yy.^2 ./ (0.25*D_t^2));
eps_c   = 10;
OF = 2.24;
p_c = 50e5;

r_cc = (sqrt(eps_c)*D_t)/2;
l_con = (r_cc - D_t/2)/tan(pi/6);
x_cnv = linspace(0,l_con,10)';
y_cnv = D_t/2 + tan(pi/3).*(l_con - x_cnv);


