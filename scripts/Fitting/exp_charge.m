function f = exp_charge(x, t)
% Exponential charge
f = x(1)*(1-exp(-t/x(2)));