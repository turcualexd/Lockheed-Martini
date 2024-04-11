function h_g = h_g_BARTZ(mu, c_p, Pr, D_t, p_c, c_star, R, sigma, A)
% Bartz model for the heat transfer coefficient

a = (0.026/D_t^0.2) * (mu^0.2*c_p/Pr^0.6) * (p_c/c_star)^0.8 * (D_t/R)^0.1 .* sigma;
A_t = pi*D_t^2 / 4;
h_g = a.*(A_t/A).^(0.9);

end