function [Bx, By, Bz, Bm] = reshape_mag(B, x, y, z)
nx = length(x);
ny = length(y);
nz = length(z);

Bx = reshape(B(:, 1), nx, ny, nz);
By = reshape(B(:, 2), nx, ny, nz);
Bz = reshape(B(:, 3), nx, ny, nz);

Bm = (Bx.^2 + By.^2 + Bz.^2).^(1/2);