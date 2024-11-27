function f = moody(Re)
if Re < 2300
    f = 64/Re;
else
    opt = optimoptions("fsolve", "Display","none");
    f = fsolve(@(x) 1/sqrt(x) + 2*log10(2.51/(Re*sqrt(x))), -4.5e-10*Re+0.051, opt);
end