function movie = animate_3out(in, f_name, x_lims, y_lims, z_lims, head_frac, radius_1, radius_2, square_off, color)
% Function for animating magnetization as a function of time.
%
% Requires arrow3d
%
% in = magnetization vectors, form (np, 3)
%
% f_name = filename. If one is not provided, you will be prompted.
%
% x_lims = x window limits [-lim +lim], default is min/max + 5% or 0
% y_lims = y window limits,[-lim +lim] default is min/max + 5% or 0
% z_lims = z window limits,[-lim +lim] default is min/max + 5% or 0
%
% These will be made square unless square_off is set to 1.
%
% head_frac = fraction of the arrow that is the head, default = 0.9
% radius_1 = radius of the arrow base, default = 0.005
% radius_2 = radius of the head base, default = 0.02
% color = color of the arrow, default is blue.
%
% Radii will be adjusted according to limits for consistency. 0.005 and
% 0.02 numbers are given for limits of -1 to 1 in all three dimensions.
%
%
% Usage: animate_3out(in, f_name, x_lims, y_lims, z_lims, head_frac, radius_1, radius_2, square_off, color);

if nargin < 3 || x_lims == 0
    x_lims = zeros(1, 2);
    [x_lims(1) x_lims(2)] = calculate_lims(in(:, 1));
end

if nargin < 4 || y_lims == 0
    y_lims = zeros(1, 2);
    [y_lims(1) y_lims(2)] = calculate_lims(in(:, 2));
end


if nargin < 5 || z_lims == 0
    z_lims = zeros(1, 2);
    [z_lims(1) z_lims(2)] = calculate_lims(in(:, 2));
end

if nargin < 6 || head_frac < 0
   head_frac = 0.9; 
end

if nargin < 7 || radius_1 <=0
    radius_1 = 0.005;
end

if nargin < 8 || radius_2 <= 0
    radius_2 = 0.02;
end

if nargin < 9
    square_off = 0;
end

if nargin < 10
    color = 'blue';
end

if nargin < 2 || ~isstr(f_name)
   [file, path] = uiputfile();
   if file == 0
       return
   end
   
   if path(end) ~= filesep
       path = strcat(path, filesep);
   end
   
   f_name = sprintf('%s%s.mat', path, file);
end


%Let's try and make this thing square - hopefully that won't be horrible.
ma = max([x_lims(2), y_lims(2), z_lims(2)]);
mi = -max([abs(x_lims(1)), abs(y_lims(1)), abs(z_lims(1))]);
    
    if 0
    %If one of the limits is zero, respect that, otherwise, max it out.
    
        if x(1) ~= 0
            x(1) = mi;
        end

        if x(2) ~= 0
            x(2) = ma;
        end

        if y(1) ~= 0
            y(1) = mi;
        end

        if y(2) ~= 0
            y(2) = ma;
        end

        if z(1) ~= 0
            z(1) = mi;
        end

        if z(2) ~= 0
            z(2) = ma;
        end
    end
    
if square_off ~= 1
    
    if ma < 1
        ma = 1;
    end
    
    if mi > -1
        mi = -1;
    end
    
    x_lims = [mi ma];
    y_lims = [mi ma];
    z_lims = [mi ma];
end

% Recalculate the radii.
mag = ma-mi;
radius_1 = mag*(radius_1/2);
radius_2 = mag*(radius_2/2);

% This is how many frames we'll have.
frames = size(in);
frames = frames(1);

h = ishold();

hold off

for i = 1:frames
   v = vec2c0(in(i, :));
   
   % Plot the axes.
   arrow3d([0, 1], [0, 0], [0, 0], 0.9, radius_1, radius_2/2, 'black');
   hold on
   arrow3d([0, -1], [0, 0], [0, 0], 0.9, radius_1, radius_2, 'black');
   arrow3d([0, 0], [0, 1], [0, 0], 0.9, radius_1, radius_2, 'black');
   arrow3d([0, 0], [0, -1], [0, 0], 0.9, radius_1, radius_2, 'black');
   arrow3d([0, 0], [0, 0], [0, 1], 0.9, radius_1, radius_2, 'black');
   arrow3d([0, 0], [0, 0], [0, -1], 0.9, radius_1, radius_2, 'black');
   
   arrow3d(v{:}, head_frac, radius_1, radius_2, color); 
  
   hold off
   
   
   xlim(x_lims);
   ylim(y_lims);
   zlim(z_lims);
     
   % Make the movie.
   m(i) = getframe; % Just gets the still frame each time.
end

save(f_name, 'm');
movie = m;

if h
    hold on
end

function [mi, ma] = calculate_lims(in)
    ma = max(in(:, 1));
    mi = min(in(:, 1));
    
    if(ma < 0)
        ma = 0;
    end
    
    if(mi > 0)
        mi = 0;
    end
    
    ma = ma+ma*0.05;
    mi = mi+mi*0.05;