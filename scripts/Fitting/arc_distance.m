function out = arc_distance(x)
s1 = length(x);

x2 = zeros(s1/2, 2);
x2(:, 1) = x(1:(s1/2));
x2(:, 2) = x((s1/2 + 1):end);

xyz = zeros(size(x2, 1), 3);
[xyz(:, 1) xyz(:, 2) xyz(:, 3)] = sph2cart(x2(:, 1), x2(:, 2), ones(size(x2, 1), 1));
s = abs(acos(xyz*xyz')).^3;
out = sum(sum(1./s(~diag(ones(1, size(s, 1))))));
% out = 1/out;