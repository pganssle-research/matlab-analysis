function grid = makegridpol(r, thet, z)
% Makes a grid out of the input vectors r, theta, z
%
% Usage:
% grid = makegridv(r, theta, z);

nr = length(r);
nt = length(thet);

if(exist('z', 'var'))
    nz = length(z);
    is3d = 1;
    
    n = nr*nt*nz;
    
    [ri, theti, zi] = ind2sub([nr, nt, nz], 1:n);
    
    grid = zeros(n, 3);
else
    is3d = 0;
    n = nx*ny;
    
    [ri, theti] = ind2sub([nr, nt], 1:n);
    
    grid = zeros(n, 2);
end
 
grid(:, 1) = arrayfun(@(i, j)r(i)*cos(thet(j)), ri, theti);
grid(:, 2) = arrayfun(@(i, j)r(i)*sin(thet(j)), ri, theti);

if(is3d)
   grid(:, 3) = arrayfun(@(i)z(i), zi); 
end