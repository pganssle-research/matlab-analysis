function out = calc_mag_sens(peak, noise, verbose, voltage, resistance) 
% Get magnetic field sensitivity for a given peak/noise, with the
% assumption that the voltage is in mVrms/sqrt(Hz), the response of the
% coil is 1.96 G/A. If resistance is left out, it is assumed to be 100k, if
% voltage is left out, it is assumed to be 16mV.

if nargin < 4
    voltage = 16; %Voltage in mVrms
end

if nargin < 5
    resistance = 100000; % Resistance in ohms.
end

snr = peak/noise;
current = voltage/resistance; % Current in mArms/sqrt(Hz)
field = current*3.92/1000; % Field in Grms/sqrt(Hz)

out = field*1e11/snr; % Sensitivity in fT/sqrt(Hz)

if verbose
    omag = floor(log10(out));
    db = 1;
    
    if omag < 3
        units = 'fT';
    elseif omag < 6
        units = 'pT';
        db = 1e3;
    elseif omag < 9
        db = 1e6;
        units = 'nT';
    elseif omag < 12
        db = 1e9;
        units = 'mT';
    end
    o = out/db;
    fprintf('Sensitivity: %f %s /sqrt(Hz)\n', o, units);
end