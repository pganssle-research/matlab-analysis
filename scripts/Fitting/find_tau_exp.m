function [tau, taue] = find_tau_exp(t, c)
% Find tau from curves
%
% Usage:
% [tau, taue] = find_tau_exp(t, c);

typical_values = [2.5, c(1)];
options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'TypicalX', typical_values, 'MaxFunEvals', 2000);
  
if(size(t, 2) > size(t, 1))
	t = t';
end

if(size(c, 2) > size(c, 1))
	c = c';
end

% Get the results, save only the exitflag which gives convergence.
warning('off'); %#ok;

[tau, ~, r, ~, ~, ~, J]  = lsqcurvefit(@exponential_fit, typical_values, t, c, [], [], options);
	
	% Standard error (95% confidence interval)
ci = nlparci(tau, r, 'jacobian', J);
taue = ci(:, 2)'-tau;
warning('on'); %#ok;
