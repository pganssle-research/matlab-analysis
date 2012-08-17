function [s, l] = wire_sensitivity(w, D, rho, W, d, T)

if(~exist('T', 'var'))
	T = 25; %STP
end

D = D*100; % Convert centimeters to meters;
W = W*100;
d = d*100; 

mu0 = 4*pi*1e-7; % In T*m/A or whatever.
T = T+273.15; % Convert to Kelvin.

Vw = pi*D*W^2;
kb = 1.38e-23;  

coeff = (8/(w*D))*sqrt(kb*T*rho);

Blf = coeff/sqrt(Vw);
Bhf = coeff/sqrt(1.8*pi*D*W*sqrt(2*rho/(w*mu0)));