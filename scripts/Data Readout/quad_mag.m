function [c_x, c_y, c_m, times, spans] = get_mag_curve(start, freq, num_windows, sr, window, data, asym, t_flag, no_offsets)
	% Function for getting the quadrature magnetization curve in an experiment
	% consisting of a train of 180x, 90y.
	%
	% Usage: [c_x, c_y, c_m, times, spans] = quad_mag(start, freq, num_windows, sr, window, data, asym, t_flag);
	%
	% start = Time point to start from (ms)
	% freq = frequency of pulses
	% num_windows = number of 'points' to acquire
	% window = how long to average between pulses (ms)
	% sr is the sampling rate
	% data = the experimental data
	% asym is an optional argument that allows for asymmetrical sampling of the medians
	% t_flag should be set high if you want the data transposed on the output, default is 1.
	%
	% out is an output matrix of size (number of indirect points)*finish
	% times is the times
	% spans is the windows sampled from

	
	% Get our data into readable variables
	data_size = size(data);
	np = data_size(1);
	nt = data_size(2);
	nd = length(data_size)-1;
	max_points = zeros(1, nd-1);
	if nargin < 6
		fprintf('Insufficient number of arguments\n')
		return
	elseif nargin > 8
		fprintf('Too many arguments!\n')
		return
	end
	
	
	if nargin < 7
		asym = 0.5;
	end
	
	if nargin < 8
		t_flag = 1;
    end
    
    if nargin <9
        no_offsets = 0;
    end
	
	for i=0:nd-2
		max_points(i+1) = data_size(i+3);
	end
	
	if mod(num_windows, 4) ~= 0
		fprintf('Number of windows must be an even multiple of 4\n')
		return
	end
	
	% Initialize the current point and number of indirect points.
	indir_points = 1;
	for i=1:nd-1
		indir_points = indir_points*max_points(i);
	end
	
	s = [num_windows/4 3 indir_points];
	out = zeros(s);
	times = zeros(2, num_windows/4);
	
	current_point = ones(1, nd-1);

	% Initialize the vectors to pull out the relevant information
	% Each point in ms is 1000/sr ms.
	% We want to start start ms in, so:
	start_point = round(start*sr/1000);
	
	% Now we need to know how often to sample, which is freq, in Hz.
	% So we have (1/sr) points/Hz, so we want to sample every freq/sr Hz.
	% We'll round to the nearest integer value obviously.
	inc_points = round(sr/freq);
	
	% Finally we need to make a matrix containing vectors of points we want.
	% Start with the window size in points:
	window_size = round(sr/1000 * window);
	
	% If the window size is bigger than the increment, that's an issue, so freak out.
	if window_size > inc_points
		fprintf('Window is larger than the period! Not allowed.\n')
		return
	end
	
	% Now we calculate the margins on either side in points
	margins = round((inc_points-window_size)*asym);
	
	% Finally initialize the matrix, which is of size num_windows*2
	% It will be the same for all dimension points, so this is just a 2 x num_windows matrix
	% of spans.
	windows = zeros(2, num_windows);
	for i = 1:num_windows
		% First the start of the span
		windows(1, i) = start_point+inc_points*(i-1)+margins;
		
		% Then the end
		windows(2, i) = start_point+inc_points*(i-1)+margins+window_size;
	end
	
	% Make a plot so that you can see where you sampled from.
	spans = zeros(np, 1);
	for i = 1:length(windows)
		for j = windows(1, i):windows(2, i)
			if mod(i, 2) == 1
				spans(j) = 1;
			else
				spans(j) = -1;
			end
		end
	end
	
	
	% Now we just go through and turn it into the output matrix
	% It goes as such: +x, -x, +y, -y, etc
	j = 0;
	while j < indir_points
		% This is a weird way to do this.
		j = j+1;
		
		avg_vec = zeros(num_windows);
		
		cp = num2cell(current_point);
		
		% Just sum them all and divide by nt later.
		for i = 1:nt
			for k = 1:num_windows
				avg_vec(k) = avg_vec(k)+mean(data(windows(1, k):windows(2, k), i, cp{:}));
			end
		end
		
		% The part where you divide by them.
		% Now you have the vector of the averages
		avg_vec = avg_vec/nt;
		
		% Current_mat is the current matrix
		% current_mat(1) is the "x" component
		% current_mat(2) is the "y" component
		current_mat = zeros(2, num_windows/4);
		for i=1:(num_windows)
			if mod(i, 4) == 1
				current_mat(1, ((i-1)/4)+1) = avg_vec(i) - avg_vec(i+1);
			elseif mod(i, 4) == 3
				current_mat(2, ((i-3)/4)+1) = avg_vec(i) - avg_vec(i+1);
			end
		end
		
		out(:, 1, cp{:}) = current_mat(1, :);
        out(:, 2, cp{:}) = current_mat(2, :);
			
		for d = 1:nd-1
			if current_point(d) < max_points(d)
				current_point(d) = current_point(d)+1;
				for k = 1:d-1
					current_point(k) = 1;
				end
				break
			end
		end
	end
	
	% Now that that's all done, we just need to give you the time series points
	% against which to plot these. Should be simple. We'll say that each point
	% is at the time between the two windows of which it's the average, and shift
	% everything over by a window's length. This means that we start at 
	% start_point + 1000/freq ms, and each point is 1000/freq ms apart.
	
	for i=1:num_windows/4
		times(1, i) = 2+(i-1)*4000/freq;
		times(2, i) = 4+(i-1)*4000/freq;
	end
	
	% Create the output vectors.
	
	s2 = size(out);
	s3 = [s2(1), s2(3)];
	c_x = zeros(s3);
	c_y = zeros(s3);
	c_m = zeros(s3);
	
    
	c_x(:, :) = out(:, 1, :);
	c_y(:, :) = out(:, 2, :);

    if ~no_offsets
        c_x_off = squeeze(c_x(:, end));
        c_y_off = squeeze(c_y(:, end));
        
        for n = 1:s3(1)
            c_m(n, :) = ((c_x(n, :) - c_x_off(n)).^2 + (c_y(n, :) - c_y_off(n)).^2).^(1/2);
        end
    else
        c_m(:, :) = (c_x(:, :).^2 + c_y(:, :).^2).^(1/2);
    end
    
    if t_flag == 1
		c_x = transpose(c_x);
		c_y = transpose(c_y);
		c_m = transpose(c_m);
	end

	