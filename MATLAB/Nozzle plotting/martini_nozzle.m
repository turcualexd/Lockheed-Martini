clear; close all; clc;

rt = sqrt(1/pi); 
eps = 300;
perc = 1e2;
[x, y, thn, the, De] = rao_nozzle(rt, eps, perc);
x_vet = linspace(0, x(end), 10);
eps_vet = eps_equispaced(rt, eps, perc, x_vet);
A = pi * y .^ 2;
%%

figure

hold on

grid on

grid minor

axis equal

plot(x, y, 'Color', 'k', 'LineWidth', 2);
plot(x_vet, eps_vet, 'o')

the
thn
lambda = (1+cos(deg2rad(the)))/2
%%
% tabn = readtable("thn_100.xlsx");
% 
% tabe = readtable("the_100.xlsx");
% 
% xn = table2array(tabn(:,1));
% 
% yn = table2array(tabn(:,2));
% 
% x = linspace(xn(1), 300, 1e3);
% y = NaN(size(x));
% 
% for k = 1:length(x)
%     if x(k) <= xn(end)
%         y(k) = interp1(xn, yn, x(k));
%     else
%         y(k) = yn(end-1) + (yn(end)-yn(end-1))/(xn(end)-xn(end-1))*(x(k)-xn(end-1));
%     end
% end
% 
% xe = table2array(tabe(:,1));
% 
% ye = table2array(tabe(:,2));
% 
% figure
% hold on
% grid minor
% 
% plot(xn,yn)
% plot(x, y)
% 
% legend('thn', 'the')
%%




