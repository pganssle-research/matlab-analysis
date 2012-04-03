function mat = gauss2d(sigma, center, s1, s2, bot, left, space1, space2)
% Generates a 2D gaussian, centered at 'center' with standard deviation
% 'sigma' on a grid of size s1xs2, where the coordinate system starts at
% [bot, left] and the interval is 'space'
%
% If 'Center' is a scalar, it will be centered on (center, center). If it
% is greater than a 2-vector, it will be centered on (center(1), center(2))
%
% Everything after s1 is optional, defaults are:
% s2 = s1
% bot = 0
% left = 0
% space1 = 1
% space2 = space1
%
% mat = gauss2d(sigma, center, s1, [s2, bot, left, space]);

if(isscalar(center))
    center = [center, center];
elseif (size(size(center)) > 2)
    center(3:end) = [];
end

if(nargin < 4)
    s2 = s1;
end

if(nargin < 5)
    bot = 0;
end

if(nargin < 6)
    left = 0;
end

if(nargin < 7)
    space1 = 1;
end

if(nargin < 8)
    space2 = space1;
end

[R, C] = ndgrid(bot:space1:(bot+(s1-1)*space1), left:space2:(left+(s2-1)*space2));

xc = center(1);
yc = center(2);
exponent = ((R-xc).^2 + (C-yc).^2./(2*sigma^2));
mat = exp(-exponent);
% mat = mat/sum(sum(mat));
