function D = calc_diffusion(T, T0, coeffs)
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

p = fliplr(coeffs);
T = T+273.15;
T0 = T0+273.15;

D = exp(polyval(p, (1./T - 1/T0)));