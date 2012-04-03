function [vec, os] = vectify(img)
% Function to vectorize a given image. This will take an image of size [n1,
% n2, ... nx, m] and reshape it to be of size prod([n1, ... nx], m). This
% is useful for dimensionality reduction techniques. Returns 'vec' and
% 'os', where 'os' is the original size of the image (for the reverse
% operation)
%
% [vec, os] = vectify(img);

s = size(img);
ns = prod(s(1:end-1));
vec = reshape(img, ns, s(end));
os = s;