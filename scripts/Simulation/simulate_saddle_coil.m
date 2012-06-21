function B = simulate_saddle_coil(h, r, mesh, anti, nt, metric)
% Following the Biot-Savart rule, simulate a saddle coil with radius r and
% height h at points mesh. Mesh should be a vector of size [n, 3] where n
% is the number of points in the mesh.
%
% Assumption is 1A, scales linearly.
%
% h and r are in inches, if metric flag is on, units are m
% Optional input cspace is the spacing (inches, metric - mm), default 0.01
%
% B = simulate_saddle_coil(h, r, mesh[, anti, nt, metric]);

histpath = 'simulate_saddle_coil_hist.mat';
truncate_at = 200;

if(exist(histpath, 'file'))
    load(histpath, 'hist_times', 'hist_nums', '-mat');
    
    if(exist('hist_times', 'var') && exist('hist_nums', 'var'))
       warning('off');
        pol = polyfit(hist_nums, hist_times, 1);
        total_time = polyval(pol, size(mesh, 1));
		  warning('on');
    end
end

if(~exist('total_time', 'var') || total_time < 0)
    total_time = 0.02*size(mesh, 1) + 0.03;
end

hours = floor(total_time/3600);
total_time = total_time-hours*3600;
minutes = floor((total_time)/60);
total_time = total_time-minutes*60;
seconds = total_time;
fprintf('Beginning calculation - expected to take %02d:%02d:%06.3f\n', hours, minutes, seconds);

if(~exist('anti', 'var'))
    anti =  0;
end

if(~exist('nt', 'var') || nt <= 0)
	nt = 12;
end

if(~exist('metric', 'var') || metric ~= 1)
    % Convert to meters.
	metric = 0;
end

start_time = cputime;

arc_angle = 2*pi/3; % 120 degrees
s1 = (pi-arc_angle)/2;
e1 = s1 + arc_angle;

s2 = s1 + pi;
e2 = e1 + pi;

top = h/2;
bot = -h/2;

% Make the curves start with top-left
c = struct('r', r, 'h', top, 'start_angle', s1, 'arc_length', arc_angle, 'nt', nt, 'direction', 1);
c = repmat(c, 1, 4);

% Bottom left
c(2).h = bot;
c(2).start_angle = s1;
c(2).direction = 0;

% Top right
c(3).start_angle = s2;

% Bottom right
c(4).h = bot;
c(4).start_angle = s2;

if(anti)
	c(3).direction = 1;
	c(4).direction = 0;
else
	c(3).direction = 0;
	c(4).direction = 1;
end
	
% Then do the lines - start with the one starting at the end_angle of the
% top left curve.
l = struct('start', [r*cos(e1), r*sin(e1), top], 'l', h, 'nt', nt, 'direction', 1);
l = repmat(l, 1, 4);

% The other line on the left loop.
l(2).start = [r*cos(s1), r*sin(s1), top];
l(2).direction = 0;

% Then the lines on the right loop.
l(3).start = [r*cos(s2), r*sin(s2), top];
l(4).start = [r*cos(e2), r*sin(e2), top];

if(anti)
	l(3).direction = 0;
	l(4).direction = 1;
else
	l(3).direction = 1;
	l(4).direction = 0;
end

current = 1;

B = simulate_any_coil(c, l, mesh, metric, current, 0);

total_time = cputime-start_time;

if(exist('hist_times', 'var') && exist('hist_nums', 'var'))
   hist_times = [total_time hist_times];
   hist_nums = [size(mesh, 1) hist_nums];
   
   if(length(hist_times) > truncate_at)
      hist_times = hist_times(1:truncate_at);
   end
   
   if(length(hist_nums) > truncate_at)
      hist_nums = hist_nums(1:truncate_at);
   end   
else
   hist_times = [total_time];
   hist_nums = [size(mesh, 1)];
end

save(histpath, 'hist_times', 'hist_nums');

hours = floor(total_time/3600);
total_time = total_time-hours*3600;
minutes = floor((total_time)/60);
total_time = total_time-minutes*60;
seconds = total_time;
fprintf('Elapsed time: %02d:%02d:%06.3f\n', hours, minutes, seconds);