o = zeros(4000, 4, 15);
for i = 0:14
    fname = sprintf('C:\\Data\\Diffusion\\HighGrad\\ConstantTime\\111010-DiffusionHighGradConstantTime%04d', i);
    [o(:, :, i+1), p, np, sr] = readout(0, fname);
end

c = get_mag(200, 10, 6, sr, 50, o, 0.75);
t = [1, 5, 10, 15, 20, 25, 30, 37.5, 40, 50, 54.545454, 60, 66.66667, 75, 100]*1e-3;
plot(t, c);
