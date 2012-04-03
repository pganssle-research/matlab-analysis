function ind = solenoid_inductance(gauge, len, radius, metric)
% Calculate the inductance of a solenoid from the gauge and length. 
%
% Length and radius should be in inches if metric is missing 
% or set to 0, otherwise it should be given in cm.
%
% Returns value in mH
%
% Usage:
% ind = solenoid_inductance(gauge, len, radius, metric)

if(~exist('metric', 'var') || metric == 0)
    len = len * 2.54;
    radius = radius * 2.54;
end

d = (exp(2.1104-0.11594*gauge))/10; % Formula for AWG gauge in mm, converted to cm
n = 2*len/d;
A = pi * radius^2; % Cross-sectional area in cm^2
mu0 = 4*pi*10^-6; % mu0 in mH/cm

ind = mu0 * n^2 * A / len; % Answer in mH