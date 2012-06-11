function B = simulate_hh_coil(h, r, mesh, anti, metric)
% Following the Biot-Savart rule, simulate a saddle coil with radius r and
% height h at points mesh. Mesh should be a vector of size [n, 3] where n
% is the number of points in the mesh.
%
% Assumption is 1A, scales linearly.
%
% h and r are in inches, if metric flag is on, units are m
% Optional input cspace is the spacing (inches, metric - mm), default 0.01
%
% B = simulate_coil(h, r, mesh[, anti, metric]);

histpath = 'simulate_hh_coil_hist.mat';
truncate_at = 200;

if(exist(histpath, 'file'))
    load(histpath, 'hist_times', 'hist_nums', '-mat');
    
    if(exist('hist_times', 'var') && exist('hist_nums', 'var'))
       
        pol = polyfit(hist_nums, hist_times, 1);
        total_time = polyval(pol, size(mesh, 1));
    end
end

current = 1;
start_time = cputime;

if(~exist('anti', 'var'))
    anti = 0;
end

% Make the curves
c = struct('r', r, 'h', h, 'start_angle', 0, 'arc_length', 2*pi, 'nt', 8, 'direction', 1);

c = repmat(c, 1, 2);

c(2).h = -h;
if(anti)
    c(2).direction = 0;
end

B = simulate_any_coil(c, [], mesh, metric, current, 0);

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


function cur = curve_xy(x0, y0, x, y, z1)
    z = repmat(z, size(x));
    
    cur = zeros(length(x), 3);
    
    x1 = x-x0;
    y1 = y-y0;
    
    d = ((x-x0).^2 + (y-y0).^2 + (z-z0).^2);
    d = repmat(d, 3, 1);
    
    cur(:, 1) = z1;
    cur(:, 2) = -z1;
    cur(:, 3) = y1*x1;
  
    cur = cur./d;

function cur = curve_theta(x0, y0, z0, R, thet, z)
ct = cos(thet);
st = sin(thet);

x = R*ct - x0;
y = R*st - y0;
z = z-z0;
z = repmat(z, size(thet));

cur = zeros(length(thet), 3);

d = (x.^2 + y.^2 + z.^2).^(3/2);

repmat(d, 3, 1);

cur(:, 1) = z.*ct;
cur(:, 2) = z.*st;
cur(:, 3) = x.*ct - y.*st;

cur = R*cur./d;

function lin = lin_int(x1, y1, z0, z)

d = (x1^2 + y1^2 + (z-z0).^2).^(-3/2);

lin = repmat(d, 1, 2);
lin(:, 1) = -lin(:, 1)*y1; % x
lin(:, 2) = lin(:, 2)*x1; % y
