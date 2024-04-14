function [x, y, thn, the, De, x3, y3] = rao_nozzle(rt, eps, perc)

switch perc

    case 60

        tabn = readtable("thn_60.xlsx");

        tabe = readtable("the_60.xlsx");

    case 70

        tabn = readtable("thn_70.xlsx");

        tabe = readtable("the_70.xlsx");

    case 80 

        tabn = readtable("thn_80.xlsx");

        tabe = readtable("the_80.xlsx");

    case 90

        tabn = readtable("thn_90.xlsx");

        tabe = readtable("the_90.xlsx");

    case 100

        tabn = readtable("thn_100.xlsx");

        tabe = readtable("the_100.xlsx");

    otherwise

        error('Only theese percentages: 60, 70, 80, 90, 100');

end

xn = table2array(tabn(:,1));

yn = table2array(tabn(:,2));

xe = table2array(tabe(:,1));

ye = table2array(tabe(:,2));

if eps <= xn(end)
    
    thn = interp1(xn, yn, eps);
else
    thn = yn(end-1) + (yn(end)-yn(end-1))/(xn(end)-xn(end-1))*(eps-xn(end-1));
end

if eps <= xe(end)
    
    the = interp1(xe, ye, eps);
else
    the = ye(end-1) + (ye(end)-ye(end-1))/(xe(end)-xe(end-1))*(eps-xe(end-1));
end

perc = perc / 100;

th1 = - 135 * pi / 180 : 0.01 * pi :- 90 * pi / 180;

x1 = 1.5 * rt * cos(th1); 

y1 = 1.5 * rt * sin(th1) + 1.5 * rt + rt;

th2 = - 90 * pi / 180: 0.01 * pi : (thn - 90) * pi / 180;

x2 = 0.382 * rt * cos(th2);

y2 = 0.382 * rt * sin(th2) + 0.382 * rt + rt;

Nx = x2(end); Ny = y2(end);

Ey = sqrt(eps) * rt; 

Ex = perc * ( (sqrt(eps) - 1) * rt / tan(pi / 12) );

m1 = tan(thn * pi / 180); m2 = tan(the * pi / 180); 

C1 = Ny - m1 * Nx; C2 = Ey - m2 * Ex; 

Qx = (C2 - C1) / (m1 - m2); Qy = (m1 * C2 - m2 * C1) / (m1 - m2);

t = 0 : 0.01 : 1; 

x3 = (1 - t) .^ 2 * Nx + 2 * (1 - t) .* t * Qx + t .^2 * Ex;

y3 = (1 - t) .^ 2 * Ny + 2 * (1 - t) .* t * Qy + t .^2 * Ey;

x = [x1 x2 x3]; y = [y1 y2 y3];

De = y(end) * 2;

end