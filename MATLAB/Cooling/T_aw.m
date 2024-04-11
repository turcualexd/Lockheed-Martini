function T_aw = T_aw(gamma, M_x, T_c, Pr)
  
   % gamma: heat capacity ratio of the gas at location x along nozzle
   % M_x  : mach number at location x along nozzle
   % T_c  : combustion chamber (stagnation temperature)
   % mu   : dynamic viscosity of the gas (at position x of nozzle)
   % c_p  : heat capcity at constant pressure at location x along nozzle
   % k    : thermal conductivy at location x
   
   r = Pr^(0.33);

   T_aw = T_c*( ( 1 + r*((gamma - 1)/2).*M_x.^2 ) ./ ( 1 + ((gamma - 1)/2).*M_x.^2 ) );
   
end

