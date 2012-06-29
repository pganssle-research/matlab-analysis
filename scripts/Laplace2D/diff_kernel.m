function y = diff_kernel(G, D, n, tau, gamma)
% Diffusion, for gradient G (G/cm) at diffusion coefficient D (1e-5 cm^2/s), 
% with n  CMPGs with interpulse spacing tau. Gamma is the gyromagnetic 
% ratio of the nucleus in Hz/Gauss.

A = -(2/3)*(2*pi*gamma)^2*n*tau.^3;
D = D*1e-5;

y = cell2mat(arrayfun(@(g)exp(A*D*(g).^2), G, ... 
	'UniformOutput', false));