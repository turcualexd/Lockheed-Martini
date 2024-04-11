clear; close all; clc;

M = 
sigma = @(x,M) 1./ ( (0.5*x.*(1 + M.^2 .* (gamma-1)/2) + 1/2).^(0.68) .* (1 + M.^2 .* (gamma - 1)/2 ).^0.12 );
gamma = 1.14;