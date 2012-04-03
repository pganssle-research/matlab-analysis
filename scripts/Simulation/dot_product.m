function o = dot_product(x, y, z, n, m)
% Generates the dot product between the two things you feed this.
x_n = squeeze(x(:, :, n));
x_m = squeeze(x(:, :, m));

y_n = squeeze(y(:, :, n));
y_m = squeeze(y(:, :, m));

z_n = squeeze(z(:, :, n));
z_m = squeeze(z(:, :, m));

o = x_n*x_m + y_n*y_m + z_n*z_m;