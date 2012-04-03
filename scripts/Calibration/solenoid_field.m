function field = solenoid_field(gauge, varargin)
% Calculate the field from a solenoid from the gauge and length. If not
% specified, current is 1A.
%
% Field  = mu0 * (N/L) * I;
% Returns value in gauss.
%
% Usage:
% field = solenoid_field(gauge, [current]);
% field = solenoid_field(gauge, [voltage, resistance]);

if(length(varargin) < 1) 
    current = 1;
elseif (length(varargin) == 1)
    current = varargin{1};
else
    current = varargin{1}/varargin{2};
end

d = (exp(2.1104-0.11594*gauge))/1000; % Formula for AWG gauge in mm, converted to meters

nd = 2/d; % Number of turns - all my solenoids are wrapped twice
mu0 = 4*pi*10^-3; % mu0 in Gauss * m / A

field = mu0*current*nd; % Permeability * current * turn density

