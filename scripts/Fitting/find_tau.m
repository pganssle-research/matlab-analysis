function [out, tau, taue, t, c] = find_tau(s, options)
% Adds a t1, t2, or D to to s.fit
%
% Options - a struct containing any of the following fields:
% .ncomp:   Numeric - number of predicted components.
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
	o = struct('ncomp', 1, 'type', 't1', 'navg', 6);
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

types = {'t1', 't2', 'diff'};
if(~isfield(o, 'type') || ~ischar(o.type) || isempty(find(strcmp(o.type, types), 1, 'first')))
	o.type = 't1';
end

if(~isfield(o, 'ncomp') || isempty(o.ncomp) || ~isnumeric(o.ncomp))
	o.ncomp = 1;
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
	i1 = find(p.ps.instrs.instr == 2, 1, 'first'); % If a malformed loop is found, grab the first loop.
	if(~isempty(i1))
		n = p.ps.instrs.data(i1);
	else
		error('Invalid number of cycles!');
	end
	
	i1 = find(bitand(2^12, p.ps.instrs.flags), 1, 'first');
	if(~isempty(i1))
		tau = p.ps.instrs.ts(i1)*1000; % Convert to 2x ms.
	else
		error('Invalid tau!');
	end
	
	if(~isfield(o, 'G_cal'))
		if(isfield(s.disp, 'G_cal') && s.disp.G_cal > 0)
			G_cal = s.disp.G_cal;
		else
			G_cal = 0.0498;
		end
	end
	
	V = G_cal*p.aovals';
	out.disp.G_cal = G_cal;
	
	ins = p.aodim(1)+1; % For now assume that it's the first thing varying.
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
perm(ins) = [];
perm = [ins, perm(:)];

ns = prod(p.maxsteps(perm(2:end)));

if(length(perm) > 1)
	c = permute(c, perm); % Put the varying dimension at the front.
else
	if(size(c, 1) < size(c, 2))
		c = c';
	end
end

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'MaxFunEvals', 2000, 'TolX', 1e-16, 'TolFun', 1e-22);
  
% Get the results, save only the exitflag which gives convergence.
if(is_diff)
	out.fit.V = V;
	t = V;
else
    if(size(t, 2) > size(t, 1))
        t = t';
    end
    
    out.fit.t = t;
end

nc = o.ncomp;
out.fit.c = c;
out.fit.cf = zeros(size(c));

if(is_t1 || is_t2)
	kern = @exponential_fit;
	typical_values = 2*ones(1, nc*2);
	typical_values(2:2:end) = c(1)/nc;
elseif(is_diff)
	kern = @(x, t)diffusion_fit(x, t, n, tau);
	typical_values = 2*ones(1, nc*2)*1e-5;
	typical_values(2:2:end) = c(1)/nc;
end

options.TypicalX = typical_values;

	
out.fit.tau = zeros(ns, nc);
out.fit.tauA = zeros(ns, nc);
out.fit.tauE = zeros(ns, nc);
out.fit.tauAE = zeros(ns, nc);

out.fit.c = c;
out.fit.cf = zeros(size(c));

out.fit.beta = zeros(ns, length(typical_values));

warning('off'); %#ok;

zm = zeros(size(typical_values));

for i = 1:ns
	[tau, ~, r, ~, ~, ~, J]  = lsqcurvefit(kern, typical_values, t, c(:, i)', zm, [], options);

	out.fit.tau(i, :) = tau(1:2:end);
	out.fit.tauA(i, :) = tau(2:2:end);
	out.fit.cf(:, i) = kern(tau, t);

	% Standard error (95% confidence interval)
	ci = nlparci(tau, r, 'jacobian', J);
	out.fit.tauE(i, :) = ci(1:2:end, 2)-tau(1);
	out.fit.tauAE(i, :) = ci(2, 2)-tau(2);

	out.fit.beta(i, :) = tau;
end


tau = out.fit.tau;
taue = out.fit.tauE;

warning('on'); %#ok