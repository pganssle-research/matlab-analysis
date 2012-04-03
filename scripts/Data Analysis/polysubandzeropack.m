o_av = squeeze(mean(o(:, :, :), 2));
t = (1:np)'/sr;
for i = 1:26
    pol(:, i) = polyfit(t(28:end), o_av(28:end, i), 3);
    o_av(:, i) = o_av(:, i) - polyval(pol(:, i), t);
    o_av(:, i) = cutandzeropack(o_av(:, i), 27);
end

s = fft(o_av);
f = 0:(sr/(np-1)):sr;