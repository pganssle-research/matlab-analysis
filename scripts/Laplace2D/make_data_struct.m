function ds = make_data_struct(x, y, z, std)
% Makes a data structure from 2 vectors and a matrix. 
% 
% Inputs:
% x:    First index vector (size n)
% y:    Second index vector (size m)
% z:    Matrix (size nxm)
% std1: Standard deviation of data in the x dimension (Default: 0)
% std2: Standard deviation of data in the y dimension (Default: std1)
%
%
% Outputs:
% ds -> Struct:
%       d.x -> x
%       d.y -> y
%       d.z -> z
%       d.std1 -> std1
%		  d.std2 -> std2
%
% Usage:
% ds = make_data_struct(x, y, z, std1, std2);

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

if(~exist('std', 'var'))
	% Estimate the standard deviations.
	k = floor(length(y)/3);
	i = 0:(k-1);
	
	std1 = sum(sum((z(:, 3*i+1) - 2*z(:, 3*i+2) + z(:, 3*i+3)).^2));
	std1 = sqrt((1/(6*k*length(x)))*std1);
	
	k = floor(length(x)/3);
	i = 0:(k-1);
	
	std2 = sum(sum((z(3*i+1, :) - 2*z(3*i+2, :) + z(3*i+3, :)).^2));
	std2 = sqrt((1/(6*k*length(y)))*std2);

	std = std1*std2;
end

ds = struct('x', x, 'y', y, 'z', z, 'std', std);