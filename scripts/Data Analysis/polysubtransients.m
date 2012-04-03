o_av = o;
t = (1:np)'/sr;
for i = 1:22
    % Polynomial subtraction
    pol(:, i) = polyfit(t(28:end), o_av(28:end, i), 3);
    o_av(:, i) = o_av(:, i) - polyval(pol(:, i), t);
    
    % Apodization
    o_av(:, i) = o_av(:, i).*exp(-t/5);    
end

nnp = 2^(ceil(log2(np)));
s = fft(o_av, nnp);
f = 0:(sr/(nnp-1)):sr;
plot(f(5600:6000), magnitude(s(5600:6000, :)))