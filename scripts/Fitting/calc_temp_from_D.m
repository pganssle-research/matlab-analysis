function T = calc_temp_from_D(D, T0, coeffs)
% Calculate the diffusion based on B and C in the equation
% ln(D(T)) = ln(D0) + sum_n(coeffs(n)*(1/T-1/T0).^2);
%
% Where n runs from 0 -> length(coeffs)
% coeffs(0) should be ln(D0)
%
%
% C defaults to 0;
% 
% Usage:
% D = calc_diffusion(T, D0, T0, B[, C]);

if(D < 1e-3)
	% Assume the units are different.
	D = D*1e5;
end

p = coeffs;
T0 = T0+273.15;

c = p(1) - log(D);
b = p(2);
a = p(3);

if(length(p) < 3 || p(3) == 0)
	T = (-c./b + 1/T0).^-1 - 273.15;
else
	T = ((-b - sqrt(b^2 - 4*a*c))./(2*a) + 1/T0).^-1 - 273.15;
end
