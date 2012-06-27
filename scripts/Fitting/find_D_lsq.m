function [D, out, De, c2, V] = find_D_lsq(struct, n, tau, G_cal, navg, temp, mag_cal)
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

s = struct;

if(~exist('G_cal', 'var') || isempty(G_cal) || ~isscalar(G_cal) || G_cal <= 0)
	if(isfield(struct, 'disp') && isfield(struct.disp, 'G_cal'))
		G_cal = s.disp.G_cal;
	else
		G_cal = 0.0498;		% Gradient calibration in G/(V*cm)
	end
end

if(~exist('navg', 'var'))
	navg = 6;				% Number of windows for the make_avg_any
end

if(~exist('temp', 'var'))
	temp = 40;
end

if(~exist('mag_cal', 'var') || isempty(mag_cal) || mag_cal <= 0)
	if(isfield(s, 'disp') && isfield(s.disp, 'mag_cal'))
		mag_cal = s.disp.mag_cal;
	else
		mag_cal = 836.52;
	end
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

out = s;
s.disp.mag_cal = mag_cal;
s.disp.G_cal = G_cal;

V = G_cal*struct.prog.aovals';

if(isfield(out, 'fit') && isfield(out.fit, 'c'))
	c2 = out.fit.c;
else
	c2 = mag_cal*squeeze(make_avg_any(struct.win.ac{:}, navg));
end

if(c2(1) < 0)
	c2 = -c2;
end

typical_values = [calc_Dw_from_temp(temp), c2(1)];

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'Typical', typical_values, 'MaxFunEvals', 2000);

warning('off'); %#ok
for i = 1:1
	[D, ~, r, ~, ~, ~, J] = lsqcurvefit(@(x, t)diffusion_fit(x, t, n, tau), typical_values, V, c2(i, :), [], [], options);

	ci = nlparci(D, r, 'jacobian', J);

	De = ci(:, 2)'-D;
	
	out.fit.D(i) = D(1);
	out.fit.DA(i) = D(2);
	
	out.fit.De(i) = De(1);
	out.fit.DeA(i) = De(2);
end
warning('on'); %#ok

out.fit.V = V;
out.fit.c = c2;
out.fit.cf = diffusion_fit(D, V, n, tau);
out.fit.D = D(1);
out.fit.DA = D(2);
out.fit.De = De(1);
out.fit.DAe = De(2);
out.fit.n = n;
out.fit.tau = tau;

D = D(1);
De = De(1);