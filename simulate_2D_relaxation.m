function out = simulate_2D_relaxation(t, G, T2, D, n, tau, noise)  
% Generates a data set simulating a 2D T2-Diffusion experiment.
%
% Inputs:
% t: Time vector (s) (T2 dimension)
% G: Gradient strength (G/cm) - Diffusion dimension
% T2: T2 relaxation constant (s) (default: 2.5s)
% D: Diffusion coefficient (cm^2/s) (default: 3.2e-5)
% n: Number of repetitions (default: 24)
% tau: Interpulse spacing in ms (default: 30)
% noise: Noise level (default: 0)
%
% Usage:
% out = simulate_2D_relaxation(t, G, T2, D, n, tau, noise);

if(~exist('t', 'var'))
	error('Must provide time vector.');
end

if(~exist('G', 'var'))
	error('Must provide gradient vector.');
end

if(~exist('T2', 'var'))
	T2 = 2.5;
end

if(~exist('D', 'var'))
	D = 3.2e-5;
end

if(~exist('n', 'var') || isempty(n) || ~isnumeric(n) || n <= 0)
	n = 24;
end

if(~exist('tau', 'var') || isempty(tau) || ~isnumeric(tau) || tau <= 0)
	tau = 30e-3;
else
	tau = tau*1e-3;
end
	
% Make sure time is a column vector and G is a row vector.
if(size(t, 2) > size(t, 1))
	t = t';
end

if(size(G, 1) > size(G, 2));
	G = G';
end

hg = 4257*2*pi;

out = exp(-t/T2)*exp(-(n*D*tau*(hg*tau*G).^2)/3);

	
	