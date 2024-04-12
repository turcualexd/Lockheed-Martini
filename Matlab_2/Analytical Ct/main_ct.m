clear, clc, close all;

syms eps pr k Ct real

eq_Ct = Ct == sqrt( 2*k^2/(k-1) * (2/(k+1)) ^ ((k+1)/(k-1)) * (1 - pr^((k-1)/k)) ) ...
                + pr*eps;
% printLatex(eq_Ct)

eq_eps = 1/eps == ((k+1)/2)^(1/(k-1)) * pr^(1/k) ...
                    * sqrt( (k+1)/(k-1) * (1 - pr^((k-1)/k)) );
% printLatex(eq_eps)

%% Solution

k = 1.24;
n = 10;
eps_vec = linspace(30, 300, n).';

pr_solve = subs(eq_eps);
Ct_solve = subs(rhs(eq_Ct));

Ct_vec = nan(n,1);
for i = 1:n
    pr_temp = vpasolve( subs(pr_solve,eps,eps_vec(i)), pr, 1e-4 );
    Ct_vec(i) = double(subs(Ct_solve, pr, pr_temp));
end

figure
plot(eps_vec, Ct_vec)