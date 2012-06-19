function D = find_D_lsq(struct, n, tau, G_cal, navg, temp)
% Feed this an mc_struct, it will calculate the diffusion coefficient from
% a least squares fit.
%
% struct:	An mc_struct containing a diffusion measurement
% G_cal:		The gradient calibration in G/(V*cm) (default: 0.0512)
% n:			Number of echos
% tau:		Interpulse spacing.
% navg:		Number of windows to average (default: 6)
% temp:		Temperature of the sample, in C (default: 40)
%
% Usage:
% D = find_D_lsq(struct, G_cal, navg, temp)

if(~exist('G_cal', 'var'))
	G_cal = 0.0512;		% Gradient calibration in G/(V*cm)
end

if(~exist('navg', 'var'))
	navg = 3;				% Number of windows for the make_avg_any
end

if(~exist('temp', 'var'))
	temp = 40;
end

if(~exist('n', 'var') || n <= 0)
	i1 = find(struct.prog.ps.instrs.instr == 2, 1, 'first'); % If a malformed loop is found, grab the first loop.
	
	if(~isempty(i1))
		n = struct.prog.ps.instrs.data(i1);
	else
		error('Invalid number of cycles!');
	end	
end

if(~exist('tau', 'var') || tau <= 0)
	i1 = find(bitand(2^12, struct.prog.ps.instrs.flags), 1, 'first');
	if(~isempty(i1))
		tau = struct.prog.ps.instrs.ts(i1)*1000; % Convert to 2x ms.
	else
		error('Invalid tau!');
	end
end

V = G_cal*struct.prog.aovals';
c2 = squeeze(make_avg_any(struct.win.ac{:}, navg));

if(c2(1) < 0)
	c2 = -c2;
end

typical_values = [calc_Dw_from_temp(temp), c2(1)];

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'Typical', typical_values, 'MaxFunEvals', 2000);

warning('off'); %#ok

D = lsqcurvefit(@(x, t)diffusion_fit(x, t, n, tau), typical_values, V, c2, [], [], options);

warning('on'); %#ok

