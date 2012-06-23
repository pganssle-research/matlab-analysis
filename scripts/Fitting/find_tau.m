function [out, tau, taue, t, c] = find_tau(s, options)
% Adds a t1, t2, or D to to s.fit
%
% Options - a struct containing any of the following fields:
% .func:		'exp' (Default)
%				'biexp'
%				'diff'
%
% .type:		't1' (Default)
%				't2'
%				'diff'
%
% .new      Boolean (Default: False)
%
% .mag_cal:	Magnetic field calibration (pT/V)
% .G_cal: Gradient calibration (G/(V*cM))
% .navg: Number of windows to average (Default: 6)
%
% Usage:
% [out, t1, t, c] = find_tau(s, options);

if(exist('options', 'var'))
	o = options;
else
	o = struct('func', 'exp', 'type', 't1', 'navg', 6);
end

if(~isfield(o, 'mag_cal') || isempty(o.mag_cal) || o.mag_cal <= 0)
	if(isfield(s, 'disp') && isfield(s.disp, 'mag_cal'))
		o.mag_cal = s.disp.mag_cal;
	else
		o.mag_cal = 836.52;
	end
end

if(~isfield(o, 'new'))
	o.new = 0;
end

funcs = {'exp', 'biexp', 'diff'};
if(~isfield(o, 'func') || ~ischar(o.func) || isempty(find(strcmp(o.func, funcs), 1, 'first')))
	o.func = 'exp';
end

types = {'t1', 't2', 'diff'};
if(~isfield(o, 'type') || ~ischar(o.type) || isempty(find(strcmp(o.type, types), 1, 'first')))
	o.type = 't1';
end

if(~isfield(o, 'navg') || o.navg <= 0)
	o.navg = 6;
end

out = s;
mag_cal = o.mag_cal;
navg = o.navg;

% Simplify function thing.
type = find(strcmp(o.type, types), 1, 'first');
is_t1 = (type == 1);
is_t2 = (type == 2);
is_diff = (type == 3);

if(~is_diff)
	func = find(strcmp(o.func, funcs), 1, 'first');
	is_exp = (func == 1);
	is_biexp = (func == 2);
else
	is_exp = 0;
	is_biexp = 0;
end

% Find the dimension we're interested in.
p = s.prog;
if(is_t1)
	ins = find(p.vtypes == 1);
	dels = [p.vdel{ins}];
	[~, i] = max(dels(end, :));
	
	ins = ins(i);
	t = p.vdel{ins}; % Time vector.
elseif(is_t2)
	p = s.prog;
	ins = find(p.ps.instrs.instr(p.vins(p.vtypes == 2)+1) == 2);
	dats = [p.vdata{ins}];
	[~, i] = max(dats(end, :));
	ins = ins(i);

	t = p.vdata{ins}; % Time vector.

	% Now convert it into time
	span = find_loop_locs(p.ps, ins);
	instrs = p.ps.instrs;
	instrs.data(span(1, 1)) = 1;

	t = t*calc_span_length(instrs, span(1, :));
elseif(is_diff)
	i1 = find(struct.prog.ps.instrs.instr == 2, 1, 'first'); % If a malformed loop is found, grab the first loop.
	if(~isempty(i1))
		n = struct.prog.ps.instrs.data(i1);
	else
		error('Invalid number of cycles!');
	end
	
	i1 = find(bitand(2^12, struct.prog.ps.instrs.flags), 1, 'first');
	if(~isempty(i1))
		tau = struct.prog.ps.instrs.ts(i1)*1000; % Convert to 2x ms.
	else
		error('Invalid tau!');
	end
	
	if(~isfield(o, 'G_cal'))
		if(isfield(s.disp, 'G_cal') && s.disp.G_cal > 0)
			G_cal = 0.0498;
		end
	end
	
	V = G_cal*struct.prog.aovals';
	out.disp.G_cal = G_cal;
else
	error('Wrong function.');
end

if(~o.new && isfield(out, 'fit') && isfield(out.fit, 'c'))
	c = out.fit.c;
else
	c = squeeze(make_avg_any([s.win.ac{:}], navg))*mag_cal;
end

if(c(1) < 0)
	c = -c;
end

out.disp.mag_cal = mag_cal;

% Get the t1 for each measurement.
perm = 1:p.nDims;
perm(p.vinsdim(ins)) = [];
perm = [p.vinsdim(ins), perm(:)];

ns = prod(p.maxsteps(perm(2:end)));

if(length(perm) > 1)
	c = permute(c, perm); % Put the varying dimension at the front.
else
	if(size(c, 1) < size(c, 2))
		c = c';
	end
end

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'MaxFunEvals', 2000);
  
% Get the results, save only the exitflag which gives convergence.
if(is_diff)
	out.fit.V = V;
else
	out.fit.t = t;
end

out.fit.c = c;
out.fit.cf = zeros(size(c));

if(is_t1 || is_t2)
	warning('off'); %#ok;
	
	if(is_exp)
		out.fit.tau = zeros(ns, 1);
		out.fit.tauA = zeros(ns, 1);
		
		out.fit.tauE = zeros(ns, 1);
		out.fit.tauAE = zeros(ns, 1);
		
		out.fit.beta = zeros(ns, 2);
		
		typical_values = [2, c(1)];
		options.TypicalX = typical_values;
		
		for i = 1:ns
			[tau, ~, r, ~, ~, ~, J]  = lsqcurvefit(@exponential_fit, typical_values, t, c(:, i), [], [], options);

			out.fit.tau(:, i) = tau(1);
			out.fit.tauA(:, i) = tau(2);
			out.fit.cf(:, i) = exponential_fit(tau, t);

			% Standard error (95% confidence interval)
			ci = nlparci(tau, r, 'jacobian', J);
			out.fit.tauE(i) = ci(1, 2)-tau(1);
			out.fit.tauAE(i) = ci(2, 2)-tau(2);
			
			out.fit.beta(i, :) = tau;
		end
	else
		out.fit.tau = zeros(ns, 2);
		out.fit.tauA = zeros(ns, 2);
		
		out.fit.tauE = zeros(ns, 2);
		out.fit.tauAE = zeros(ns, 2);
		
		out.fit.beta = zeros(ns, 4);
	
		typical_values = [1, c(1)/2, 3, c(1)/2];
		
		options.TypicalX = typical_values;
		
		for i = 1:ns
			[tau, ~, r, ~, ~, ~, J]  = lsqcurvefit(@bi_exp_fit, typical_values, t, c(:, i), [], [], options);

			out.fit.tau(i, :) = tau([1 3]);
			out.fit.tauA(i, :) = tau([2 4]);
			out.fit.cf(:, i) = bi_exp_fit(tau, t);

			% Standard error (95% confidence interval)
			ci = nlparci(tau, r, 'jacobian', J);
			out.fit.tauE(i, :) = ci([1 3], 2)'-tau([1 3]);
			out.fit.tauAE(i, :) = ci([2 4], 2)'-tau([2 4]);
			
			out.fit.beta(i, :) = tau;
		end
	end
	
	tau = out.fit.tau;
	taue = out.fit.tauE;
	
	warning('on'); %#ok
else
	
end