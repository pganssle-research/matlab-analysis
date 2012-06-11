function B = simulate_line(h, r, mesh, anti, metric);
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

if(~exist('anti', 'var'))
    anti =  0;
end

if(~exist('metric', 'var') || metric ~= 1)
    % Convert to meters.
    h = h * 0.0254;
    r = r * 0.0254;
    mesh = mesh * 0.0254;
end

arc_angle = 2*pi/3; % 120 degrees
start_angle1 = (pi-arc_angle)/2;
end_angle1 = start_angle1 + arc_angle;

top = h/2;
bot = -h/2;

% Now we'll do the line parts.

llbt = zeros(size(mesh)); % Line on the left, bottom to top
lltb = zeros(size(mesh)); % Line on the left, top to bottom;
lrbt = zeros(size(mesh)); 
lrtb = zeros(size(mesh));

% Line on the left, top to bottom - at end_angle1
lltb(:, 1:2) = cell2mat(arrayfun(@(x1, y1, z0)quadv(@(z)lin_int(x1, y1, z0, z), top, bot), r*cos(end_angle1)-mesh(:, 1), ...
    r*sin(end_angle1)-mesh(:, 2), mesh(:, 3), 'UniformOutput', false));

B = lltb;

mu = 1e-3; % Vacuum permeability in meters * Gauss / A.
B = B*mu;


function cur = curve_theta(x0, y0, z0, R, thet, z)
x = R*cos(thet);
y = R*sin(thet);
z = repmat(z, size(thet));

cur = zeros(length(thet), 2);
d = (((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(3/2));
repmat(d, 2, 1);

%cur(1, :) = sqrt((x-x0).^2 + (y-y0).^2)./(((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(3/2));
%cur(2, :) = (z-z0)./(((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(3/2));

cur(:, 1) = sqrt((x-x0).^2 + (y-y0).^2);
cur(:, 2) = (z-z0);

cur = cur./d;

function lin = lin_int(x1, y1, z0, z)

a = x1^2 + y1^2;
d = (a + (z-z0).^2).^(3/2);

lin = repmat(d, 1, 2);
lin(:, 1) = lin(:, 1)*x1; % x
lin(:, 2) = lin(:, 2)*y1; % y

