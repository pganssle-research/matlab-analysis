function tau = find_tau_biexp(curves, indir_times, num_measurements, start_index, end_index)
% Function for getting the t1 from these curves
%
% Usage:
% [t1, t1_std, t1s] = find_t1(curves, indir_times, num_measurements, start_index, end_index);
%
% Only the first two arguments are mandatory.

% Extract the defaults
s = size(curves);
st = size(indir_times);

if(st(2) > st(1))
	indir_times = indir_times';
	if(find(st(2) == s) == 2)
		curves = curves';
		
		s = size(curves);
	end
end

if nargin < 3 || num_measurements == 0
    num_measurements = s(2);
end

if  nargin < 4 || start_index == 0
    start_index = 1;
end

if nargin < 5 || end_index == 0
    end_index = s(1);
end

if nargin < 2
    fprintf('Insufficient number of arguments\n')
    return
end

if nargin > 5
    fprintf('Too many arguments\n')
    return
end

% Get the t1 for each measurement.
fit = zeros(num_measurements, 3);

typical_values = [1.61, 43.5/2, 3.39, 76.6/2];
options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'Typical', typical_values, 'MaxFunEvals', 2000);

for i = 1:num_measurements
    v = curves(:, i);
    v = v(start_index:end_index);
    
    indir_times = indir_times(start_index:end_index);
    
    % Get the results, save only the exitflag which gives convergence.
    warning('off'); %#ok;
    [x, ~, ~, exitflag] = lsqcurvefit(@bi_exp_fit, typical_values, indir_times, v, [], [], options);
    warning('on'); %#ok
    fit(i, 1:2) = x(1:2);
    fit(i, 3) = exitflag;
end

tau = x;