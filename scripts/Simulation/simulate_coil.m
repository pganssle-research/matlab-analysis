function B = simulate_coil(h, r, mesh, anti, metric)
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

histpath = 'simulate_coil_hist.mat';
truncate_at = 200;

if(exist(histpath, 'file'))
    load(histpath, 'hist_times', 'hist_nums', '-mat');
    
    if(exist('hist_times', 'var') && exist('hist_nums', 'var'))
       
        pol = polyfit(hist_nums, hist_times, 1);
        total_time = polyval(pol, size(mesh, 1));
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

start_time = cputime;

if(~exist('metric', 'var') || metric ~= 1)
    % Convert to meters.
    h = h * 0.0254;
    r = r * 0.0254;
    mesh = mesh * 0.0254;
end

arc_angle = 2*pi/3; % 120 degrees
start_angle1 = (pi-arc_angle)/2;

end_angle1 = start_angle1 + arc_angle;

if(anti)
    start_angle2 = start_angle1 + pi;
    end_angle2 = end_angle1+pi;
else
    start_angle2 = end_angle1 + pi;
    end_angle2 = start_angle1 + pi;
end

start_angle = start_angle1;
end_angle = end_angle1;

top = h/2;
bot = -h/2;

%top_left = zeros(size(mesh));
%bottom_left = zeros(size(mesh));
%top_right = zeros(size(mesh));
%bottom_right = zeros(size(mesh));

top_left = cell2mat(arrayfun(@(x0, y0, z0)quadv(@(thet)curve_theta(x0, y0, z0, r, thet, top), start_angle, end_angle), ...
                        mesh(:, 1), mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

bottom_left = cell2mat(arrayfun(@(x0, y0, z0)quadv(@(y)curve_theta(x0, y0, z0, r, y, bot), end_angle, start_angle), ...
                        mesh(:, 1), mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

%[bottom_left(:, 1), bottom_left(:, 2), bottom_left(:, 3)] = pol2cart(bottom_left(:, 1), bottom_left(:, 2), bottom_left(:, 3));
   
% Reverse the angles for the next coil.
start_angle = start_angle2;
end_angle = end_angle2;

top_right = cell2mat(arrayfun(@(x0, y0, z0)quadv(@(y)curve_theta(x0, y0, z0, r, y, top), start_angle, end_angle), ...
                        mesh(:, 1), mesh(:, 2), mesh(:, 3), 'UniformOutput', false));
%[top_right(:, 1), top_right(:, 2), top_right(:, 3)] = pol2cart(top_right(:, 1), top_right(:, 2), top_right(:, 3));

bottom_right = cell2mat(arrayfun(@(x0, y0, z0)quadv(@(y)curve_theta(x0, y0, z0, r, y, bot), end_angle, start_angle), ...
                        mesh(:, 1), mesh(:, 2), mesh(:, 3), 'UniformOutput', false));
%[bottom_right(:, 1), bottom_right(:, 2), bottom_right(:, 3)] = pol2cart(bottom_right(:, 1), bottom_right(:, 2), bottom_right(:, 3));

% Now we'll do the line parts.

llbt = zeros(size(mesh)); % Line on the left, bottom to top
lltb = zeros(size(mesh)); % Line on the left, top to bottom;
lrbt = zeros(size(mesh)); 
lrtb = zeros(size(mesh));

% Line on the left, top to bottom - at end_angle1
lltb(:, 1:2) = cell2mat(arrayfun(@(x1, y1, z0)quadv(@(z)lin_int(x1, y1, z0, z), top, bot), r*cos(end_angle1)-mesh(:, 1), ...
    r*sin(end_angle1)-mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

% Line on the left, bottom to top - at start_angle1
llbt(:, 1:2) = cell2mat(arrayfun(@(x1, y1, z0)quadv(@(z)lin_int(x1, y1, z0, z), bot, top), r*cos(start_angle1)-mesh(:, 1), ...
    r*sin(start_angle1)-mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

% Line on the left, top to bottom - at end_angle2
lrtb(:, 1:2) = cell2mat(arrayfun(@(x1, y1, z0)quadv(@(z)lin_int(x1, y1, z0, z), top, bot), r*cos(end_angle2)-mesh(:, 1), ...
    r*sin(end_angle2)-mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

% Line on the left, bottom to top - at start_angle2
lrbt(:, 1:2) = cell2mat(arrayfun(@(x1, y1, z0)quadv(@(z)lin_int(x1, y1, z0, z), bot, top), r*cos(start_angle2)-mesh(:, 1), ...
    r*sin(start_angle2)-mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

B = bottom_left+top_right+bottom_right+top_left;
B = B + lltb;
B = B + llbt;
B = B + lrtb;
B = B + lrbt;

mu = 1e-3; % Vacuum permeability in meters * Gauss / A.
B = B*mu;

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

d = (x.^2 + y.^2 + z.^2).^(-3/2);

cur = R*repmat(d, length(thet), 3);

cur(:, 1) = cur(:, 1).*z.*ct;
cur(:, 2) = cur(:, 2).*z.*st;
cur(:, 3) = -cur(:, 3).*(x.*ct + y.*st);

function lin = lin_int(x1, y1, z0, z)

d = (x1^2 + y1^2 + (z-z0).^2).^(-3/2);

lin = repmat(d, 1, 2);
lin(:, 1) = -lin(:, 1)*y1; % x
lin(:, 2) = lin(:, 2)*x1; % y
