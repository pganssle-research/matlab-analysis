function grid = makegridv(x, y, z)
% Makes a grid out of the input vectors x, y and z
%
% Usage:
% grid = makegridv(x, y, z);

nx = length(x);
ny = length(y);

if(exist('z', 'var'))
    nz = length(z);
    is3d = 1;
    
    n = nx*ny*nz;
    
    [xi, yi, zi] = ind2sub([nx, ny, nz], 1:n);
    
    grid = zeros(n, 3);
else
    is3d = 0;
    n = nx*ny;
    
    [xi, yi] = ind2sub([nx, ny], 1:n);
    
    grid = zeros(n, 2);
end

grid(:, 1) = arrayfun(@(i)x(i), xi);
grid(:, 2) = arrayfun(@(i)y(i), yi);

if(is3d)
   grid(:, 3) = arrayfun(@(i)z(i), zi); 
end