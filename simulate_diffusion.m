function out = simulate_diffusion(G, n, tau, D, mag)
% Simulates a diffusion measurement for given values of the gradient G.
%
% G = gradient strength in Gauss/cm
% n = Number of repetitions
% tau = interpulse spacing, in ms
% D = diffusion coefficient
% mag = initial magnitude (default: 1)
%
% Usage:
% simulate_diffusion(G, n, tau, D[, mag])

hg = 4257*(2*pi);

if(~exist('mag', 'var'))
   mag = 1; 
end

tau = tau*1e-3;

out = mag*exp(-(((hg*G*tau).^2)*D*n*tau/3));
