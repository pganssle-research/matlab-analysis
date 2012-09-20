function [out, t, c] = load_curves(s, options)
% Retrieves the curves from the data in either a direct or indirect
% dimension.
%
% Options:
% .type    't1'
%			  't2'
%          'diff'
%          'field'
%			  'direct'
%
% .mag_cal: Magnetic field calibration
% .G_cal: Gradient calibration
% .Z_cal: Z field voltage calibration
% .n_avg: Number of windows to average
% .off: Offset
%
% [out, t, c] = load_curves(s, options);

if(~exist('options', 'var'))
	o = options;
else
	error('Must submit an options structure.');
end

if(~is_field(o, 'type') || isempty(o.type))
	error('Must choose at least one valid type!')
end

if(~is_field(o, 'mag_cal' || isempty(o.mag_cal) || o.mag_cal <= 0))
	if(isfield(s, 'disp') && isfield(s.disp, 'mag_cal'))
		o.mag_cal = s.disp.mag_cal;
	else
		o.mag_cal = 1;
	end
end

types = {'t1', 't2', 'diff', 'field', 'direct'};