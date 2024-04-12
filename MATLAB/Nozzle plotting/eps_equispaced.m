function eps_vet = eps_equispaced(rt, eps, perc, x)

eps_vet = NaN(size(x));

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

x_thn = 0.382 * rt * cos((thn - 90) * pi / 180);

Nx = x_thn;

Ny = 0.382 * rt * sin((thn - 90) * pi / 180) + 0.382 * rt + rt;

Ey = sqrt(eps) * rt;

Ex = perc * ( (sqrt(eps) - 1) * rt / tan(pi / 12) );

m1 = tan(thn * pi / 180); m2 = tan(the * pi / 180);

C1 = Ny - m1 * Nx;

C2 = Ey - m2 * Ex; 

Qx = (C2 - C1) / (m1 - m2);

Qy = (m1 * C2 - m2 * C1) / (m1 - m2);

a = Nx -2 * Qx + Ex; 

b = 2 * (Qx - Nx);

for k = 1 : length(x)

    if x(k) <= x_thn

        th = acos(x(k)./(0.382*rt));

        eps_vet(k) = 0.382 * sin(th) + 0.382*rt + rt;

    else

        c = Nx - x(k);

        t = (-b+sqrt(b^2-4*a*c))/(2*a);

        if t < 0 || t > 1

            err('t non valido!')

        end

        eps_vet(k) = (1 - t) .^ 2 * Ny + 2 * (1 - t) .* t * Qy + t .^2 * Ey;

    end

end




