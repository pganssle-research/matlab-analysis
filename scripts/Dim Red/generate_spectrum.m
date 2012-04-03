function [spec, p] = generate_spectrum(peaks, nH, np, lw)

if(~exist('lw', 'var'))
    lw = 0.002;
end

p = linspace(-12, 12, np);


spec = zeros(1, np);

for i = 1:length(peaks)
    spec = spec + lorentz(p, peaks(i), lw)*nH(i);
end


function l = lorentz(p, p0, lw)
l = (1/pi)*(lw/2)./((p-p0).^2 + (lw/2).^2);


