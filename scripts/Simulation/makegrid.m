function grid = makegrid(xs, xe, nx, ys, ye, ny, zs, ze, nz)
% Create a grid of size (nx * ny * nz) x 3 which contains all points
% linearly spaced from xs->xe, ys->ye, zs->ze in nx, ny, and nz points;
%
% Usage:
% grid = makegrid(xs, xe, nx, ys, ye, ny, zs, ze, nz);

if(nargin < 7)
    % It's 2D
    s = [nx*ny, 2];
    n = nx*ny;
    grid = zeros(n, 2);
    [xi, yi] = ind2sub([nx, ny], 1:n);
    is3d = 0;
else
    n = nx*ny*nz;
    grid = zeros(n, 3);
    [xi, yi, zi] = ind2sub([nx, ny, nz], 1:n);
    is3d = 1;
end

x = linspace(xs, xe, nx)';
grid(:, 1) = arrayfun(@(i)x(i), xi);


y = linspace(ys, ye, ny)';
grid(:, 2) = arrayfun(@(i)y(i), yi);

if(is3d)
    z = linspace(zs, ze, nz)';
    grid(:, 3) = arrayfun(@(i)z(i), zi);
end