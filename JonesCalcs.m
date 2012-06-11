function [t, out, f, s] = JonesCalcs(freq, np, sr)
% Signal for a precessing signal at frequency "freq"
%
% Usage:
% out = JonesCalcs(freq, np, t);

% Start with horizontally polarized stuff.
pol = [1; 0];

pol = repmat(pol, 1, np);
t = linspace(0, np/sr, np);
thet = 0;
pemfreq = 50000;
bet = pi/2;

for i = 1:np
    % Rotate with half-waveplate
    pol(:, i) = birefringence(pi, pi/5)*pol(:,i);
    
    % Rotate with the cell.
    %pol(:, i) = pol(:, i)+0.01*pol(:, i)*sin(2*pi*60);
    
    %pol(:, i) = rot((pi/32)*cos(2*pi*cos(2*pi*freq*t(i))))*pol(:, i);
    
    % Quarter waveplate at angle thet.
    %pol(:, i) = birefringence(-pi/2, thet)*pol(:, i);
    
    % PEM, frequency is 50k.
    pol(:, i) = birefringence(-bet*cos(2*pi*pemfreq*t(i)), pi/4)*pol(:, i);
    
    % Linear polarizer
    pol(:, i) = lp(pi/4)*pol(:, i);
end

T = 4096;

out = zeros(1, np-T);
for i = (T+1):np
    out(i-T) = (1/T)*sum(-cos(2*pi*pemfreq*t((i-T):i)).*squeeze(pol(1, (i-T):i)));
end

t = t((T+1):end);

% out = sum(pol);

npfft = 2^(floor(log2(np))+1);
s = fftx(out, npfft);
f = linspace(0, sr/2, npfft/2);
s = s(1:(npfft/2));

function out = lp(theta)
out = rot(theta)*[1, 0; 0 0];

function out = rot(theta)
out = [cos(theta), -sin(theta); sin(theta), cos(theta)];

