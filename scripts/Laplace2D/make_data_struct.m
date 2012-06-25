function ds = make_data_struct(x, y, z)
% Makes a data structure from 2 vectors and a matrix. 
% 
% Inputs:
% x:    First index vector (size n)
% y:    Second index vector (size m)
% z:    Matrix (size nxm)
%
% Outputs:
% ds -> Struct:
%       d.x -> x
%       d.y -> y
%       d.z -> z
%
% Usage:
% ds = make_data_struct(x, y, z);

if(~isvector(x) || ~isvector(y))
   error('X and Y must be vectors.'); 
end

if(size(x, 1) < size(x, 2))
   x = x';
end

if(size(y, 1) < size(y, 2))
    y = y';
end

sz = size(z);
if(~(length(x) == sz(1)))
    if(length(x) == sz(2))
        z = z';
        sz = size(z);
    else
       error('Z must be of size nxm'); 
    end
end

if(~length(y) == sz(2))
   error('Z must be of size nxm'); 
end

ds = struct('x', x, 'y', y, 'z', z);