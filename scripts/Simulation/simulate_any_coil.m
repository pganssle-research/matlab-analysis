function B = simulate_any_coil(curves, lines, mesh, metric, current, verbose)
% Following the Biot-Savart rule, simulate a saddle coil with radius r and
% height h at points mesh. Mesh should be a vector of size [n, 3] where n
% is the number of points in the mesh.
%
% Create any shape from curves and lines, which are structs.
% Here is how to create a curve:
% c.r               => Radius
% c.start_angle     => Beginning angle (0, 90, wherever it starts)
% c.arc_length      => How long is the arc
% c.h               => Height
% c.nt              => Number of turns (default 1)
% c.direction       => Direction of current flow. If it evaluates to true,
%                      it goes forward (default forward)
%
% Here is how to create a line:
% l.start           => 3-vector, [x, y, z], the start of the line
% l.l               => Length (can be negative or positive but not 0)
% l.nt              => Number of turns (default 1)
% l.direction       => Direction of current fowlo, if it evaluates to true,
%                      it goes forward (default forward)
%
% Default options for metric and current are 1 and 1A respectively.
%
% B = simulate_coil(curves, lines, mesh[, metric, current]);

histpath = 'simulate_any_coil_hist.mat';
truncate_at = 200;
num_lines = length(lines);
num_curves = length(curves);

if(~exist('verbose', 'var'))
   verbose = 0; 
end

if(exist(histpath, 'file'))
    load(histpath, '-mat');
    
    if(verbose && exist('hbtime', 'var') && exist('hbn', 'var') ...
        && exist('hltime', 'var') && exist('hln', 'var') ...
        && exist('hctime', 'var') && exist('hcn', 'var'))
                      
        polb = polyfit(hbn, hbtime, 1); %#ok
        poll = polyfit(hln, hltime, 1); %#ok
        polc = polyfit(hcn, hctime, 1); %#ok
        
        total_time = polyval(poll, num_lines*size(mesh, 1)) + polyval(polc, num_curves*size(mesh, 1)) + polyval(polb, size(mesh, 1));
    end
end

if(verbose) 
    if(~exist('total_time', 'var') || total_time < 0)
        total_time = 0.02*size(mesh, 1) + 0.03;
    end

    hours = floor(total_time/3600);
    total_time = total_time-hours*3600;
    minutes = floor((total_time)/60);
    total_time = total_time-minutes*60;
    seconds = total_time;
    fprintf('Beginning calculation - expected to take %02d:%02d:%06.3f\n', hours, minutes, seconds);
end

if(~exist('current', 'var'))
    current = 1.0;
end

start_time = cputime;

if(~exist('metric', 'var') || metric ~= 1)
    mesh = mesh * 25.4;
    
    metric = 0;
end  
    
for i = 1:length(lines)
    l = lines(i);
    if(~isfield(l, 'start') || length(l.start) < 3)
        fprintf('Line %d is missing a well-formed start point!', i);
        B = 0;
        return;   
    end
    
   if(~isfield(l, 'l') || l.l == 0)
        fprintf('Line %d is missing a well-formed length.', i);
        B = 0;
        return;   
    end
    
    if(~isfield(l, 'nt'))
        lines(i).nt = 1;
    end
    
    if(~isfield(l, 'direction'))
        lines(i).direction = 1; 
    end
    
    if(~metric)
        lines(i).l = l.l* 25.4;
        lines(i).start = l.start*25.4;
    end
end

for i = 1:length(curves)
    c = curves(i);
    if(~isfield(c, 'r') || c.r <= 0)
        fprintf('Curve %d is missing a well-formed radius!', i);
        B = 0;
        return;
    end
   
    if(~isfield(c, 'start_angle'))
        fprintf('Curve %d is missing a start angle.', i);
        B = 0;
        return;
    end
   
    if(~isfield(c, 'arc_length') || c.arc_length <= 0 || c.arc_length > 2*pi)
        fprintf('Curve %d is missing a well-formed arc length.', i);
        B = 0;
        return;
    end

    if(~isfield(c, 'h'))
        fprintf('Curve %d is missing a height', i);
      
        B = 0;
        return;
    end
   
    if(~isfield(c, 'nt'))
        curves(i).nt = 1;
    end
    
    if(~isfield(c, 'direction'))
        curves(i).direction = 1; 
    end
    
    if(~metric)
        curves(i).r = c.r * 25.4;
        curves(i).h = c.h * 25.4;
    end
end

B = zeros(size(mesh));

ctime = cputime;
for i = 1:length(curves)
    c = curves(i);
	 
    if(c.direction)
        start_angle = c.start_angle;
        end_angle = c.start_angle+c.arc_length;
    else
        start_angle = c.start_angle+c.arc_length;
        end_angle = c.start_angle;
    end
    
    a = 3;
    
    mag_buff = c.nt*cell2mat(arrayfun(@(x0, y0, z0)quadv(@(thet)curve_theta(x0, y0, z0, c.r, thet, c.h), start_angle, end_angle), ...
                        mesh(:, 1), mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

    B = B+mag_buff;
end
ctime = cputime-ctime;

mag_buff = zeros(size(mesh));

% Now we'll do the line parts.
ltime = cputime;
for i = 1:length(lines)
	l = lines(i);
    if(l.direction)
       top = l.start(3);
       bot = l.start(3)-l.l;
    else
       top = l.start(3)-l.l;
       bot = l.start(3);
    end
        
    mag_buff(:, 1:2) = l.nt*cell2mat(arrayfun(@(x1, y1, z0)quadv(@(z)lin_int(x1, y1, z0, z), top, bot), l.start(1)-mesh(:, 1), ...
        l.start(2)-mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

    B = B + mag_buff;
end
ltime = cputime - ltime;

mu = 1; % Vacuum permeability in mm * Gauss / A.
B = B*mu*current;

total_time = cputime-start_time;
btime = total_time-(ltime+ctime);

if(verbose && exist('hbtime', 'var') && exist('hbn', 'var') ...
        && exist('hltime', 'var') && exist('hln', 'var') ...
        && exist('hctime', 'var') && exist('hcn', 'var'))
   
    hbtime = [btime hbtime];
    hbn = [size(mesh, 1) hbn];

    hltime = [ltime hltime];
    hln = [num_lines*size(mesh, 1) hln];
    
    hctime = [ctime hltime];
    hcn = [num_curves*size(mesh, 1) hcn];
    
    if(length(hbtime) > truncate_at)
        hbtime = hbtime(1:truncate_at); %#ok
    end
    
    if(length(hbn) > truncate_at)
        hbn = hbn(1:truncate_at); %#ok
    end
    
    if(length(hctime) > truncate_at)
        hctime = hctime(1:truncate_at); %#ok
    end
    
    if(length(hcn) > truncate_at)
        hcn = hcn(1:truncate_at); %#ok
    end
    
    if(length(hltime) > truncate_at)
        hltime = hltime(1:truncate_at);
    end
    
    if(length(hln) > truncate_at)
        hln = hln(1:truncate_at);
    end
    
else
    hctime = ctime;
    hltime = ltime;
    hbtime = btime;
    
    hbn = size(mesh, 1);
    hln = hbn*num_lines;
    hcn = hbn*num_curves;
end

save(histpath, 'hbn', 'hbtime', 'hcn', 'hctime', 'hln', 'hltime');

if(verbose)
    hours = floor(total_time/3600);
    total_time = total_time-hours*3600;
    minutes = floor((total_time)/60);
    total_time = total_time-minutes*60;
    seconds = total_time;
    fprintf('Elapsed time: %02d:%02d:%06.3f\n', hours, minutes, seconds);
end

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

x = R*ct-x0;
y = R*st-y0;
z = z-z0;
z = repmat(z, size(thet));

d = (x.^2 + y.^2 + z.^2).^(-3/2);

cur = R*repmat(d, length(thet), 3);

cur(:, 1) = cur(:, 1).*(z.*ct);
cur(:, 2) = cur(:, 2).*(z.*st);
cur(:, 3) = -cur(:, 3).*(x.*ct + y.*st);

function lin = lin_int(x1, y1, z0, z)

d = (x1^2 + y1^2 + (z-z0).^2).^(-3/2);

lin = repmat(d, 1, 2);
lin(:, 1) = -lin(:, 1)*y1; % x
lin(:, 2) = lin(:, 2)*x1; % y
