function D = find_D_lsq(struct, temp)
% Feed this an mc_struct, it will calculate the diffusion coefficient from
% a least squares fit.
%
% Usage:
% D = find_D_lsq(struct)

G_cal = 0.0512;		% Gradient calibration in G/(V*cm)
navg = 6;				% Number of windows for the make_avg_any

if(~exist('temp', 'var'))
	temp = 40;
end

V = 0.0512*struct.prog.aovals';
c2 = squeeze(make_avg_any(out.win.ac{:}, navg));

if(c2(1) < 0)
	c2 = -c2;
end

typical_values = [calc_Dw_from_temp(temp), c2(1)];

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'Typical', typical_values, 'MaxFunEvals', 2000);

warning('off'); %#ok

D = lsqcurvefit(@diffusion_fit, typical_values, V, c2, [], [], options);

warning('on');

