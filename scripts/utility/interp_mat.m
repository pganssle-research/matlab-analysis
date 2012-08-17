function m12 = interp_mat(m1, m2, dim)
% Interpolate two matrices along a given dimension. They must have the same
% size in the other dimensions.
%
% Default value for dim is 1;
%
% m12 = interp_mat(m1, m2[, dim]);

if(~exist('dim', 'var'))
	dim = 1;
end

m1l = size(m1, dim);
m2l = size(m2, dim);

if(m1l > m2l)
	m1(m2l:end, :) = [];
elseif(m2l > m1l)
	m2(m1l:end, :) = [];
end

if(size(m1) ~= size(m2))
	error('Incommensurate sizes');
end

s = size(m1);
s(dim) = s(dim)*2;

m12 = zeros(s);

perm = 1:length(s);

perm([1, dim]) = [dim, 1];

m12 = permute(m12, perm);
m1 = permute(m1, perm);
m2 = permute(m2, perm);

m12(1:2:end, :) = m1(:, :);
m12(2:2:end, :) = m2(:, :);

m12 = permute(m12, perm);