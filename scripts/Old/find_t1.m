function [t1, t1_std, t1s] = find_t1(curves, indir_times, num_measurements, start_index, end_index)
	% Function for getting the t1 from these curves
	
	% Extract the defaults
	s = size(curves);
	if nargin < 3 || num_measurements == 0
		num_measurements = s(1);
	end
	
	if  nargin < 4 || start_index == 0
		start_index = 1;
	end
	
	if nargin < 5 || end_index == 0
		end_index = s(2);
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
	fit = zeros(num_measurements, 1);
	for i = 1:num_measurements
		v = curves(i, :);
	
		% Subtract off the minimum
		% This may be a mistake.
		v = v-min(v);
		
		v = transpose(v(start_index:end_index));
		indir_times = indir_times(start_index:end_index);
		p = polyfit(indir_times, log(v), 1);
		fit(i) = 1/p(1);
    end

    t1s = fit(fit ~= 0 & real(fit) == imag(fit) & ~isinf(fit) & ~isnan(fit));
    	
	t1s = transpose(t1s);
	if length(t1s) < 1
		fprintf('No real answers found!\n')
		return
	end
	
	t1 = mean(t1s);
	t1_std = std(t1s);