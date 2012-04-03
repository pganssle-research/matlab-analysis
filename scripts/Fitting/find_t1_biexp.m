function [t1, t1_std, t1s] = find_t1_biexp(curves, indir_times, num_measurements, start_index, end_index)
	% Function for getting the t1 from these curves
	
	% Extract the defaults
	s = size(curves);
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
	fit = zeros(num_measurements, 6);
	o = optimset('Display', 'off');
    for i = 1:num_measurements
		v = curves(:, i);
        v = v(start_index:end_index);
		
		indir_times = indir_times(start_index:end_index);
        
        % Get the results, save only the exitflag which gives convergence.
		[x resnorm resid exitflag] = lsqcurvefit(@bi_exp_fit, [3.0, 0.5, 1.0, 0.5, 0.1], indir_times, v, -[10 10 10 10 10], [10 10 10 10 10], o);   
		fit(i, 1:5) = x(1:5);
        fit(i, 6) = exitflag;
	end
	
	t1s = [];
	for i=1:num_measurements
		if fit(i, 6) >= 0;
			t1s = [t1s ; fit(i, 1:5)];
		end
	end
	
	if length(t1s) < 1
		fprintf('No real answers found!\n')
		return
	end
	
	t1 = [mean(t1s(:, 1)), mean(t1s(:, 2)), mean(t1s(:, 3)), mean(t1s(:, 4)), mean(t1s(:, 5))];
	t1_std = [std(t1s(:, 1)) std(t1s(:, 2)) std(t1s(:, 3)) std(t1s(:, 4)) std(t1s(:, 5))];