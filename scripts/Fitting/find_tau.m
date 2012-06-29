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
	types = p.vtypes(p.vtypes ~= 0);
	
	dim = p.vinsdim(types == 1);
	
	dels = [p.vdel{dim}];
	[~, i] = max(dels(end, :));
	
	dim = dim(i);
	t = [p.vdel{dim}]; % Time vector.
	
	name = 't1';
	namea = 't1a';
	
	namee = 't1e';
	nameae = 't1ae';
	
	nameb = 't1b';
elseif(is_t2)
	p = s.prog;
	
	types = p.vtypes(p.vtypes ~= 0);
		
	ins = find(p.ps.instrs.instr(p.vins(types == 2)+1) == 2);
	dim = p.vinsdim(ins);	
	
	dats = [p.vdata{dim}];
	[~, i] = max(dats(end, :));
	dim = dim(i);
	
	t = dats; % Time vector.

	% Now convert it into time
	span = find_loop_locs(p.ps, p.vins(ins));
	instrs = p.ps.instrs;
	instrs.data(span(1, 1)) = 1;

	t = t*calc_span_length(instrs, span(1, :));
	
	name = 't2';
	namea = 't2a';
	
	namee = 't2e';
	nameae = 't2ae';
	
	nameb = 't2b';
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
	
	dim = find(p.vtypes == 0, 1, 'first');
	
	name = 'D';
	namea = 'Da';
	
	namee = 'De';
	nameae = 'Dae';
	
	nameb = 'Db';
else
	error('Wrong function.');
end

if(~o.new && isfield(out, 'fit') && (isfield(out.fit, 'co') || isfield(out.fit, 'c')))
	if(isfield(out.fit, 'co'))
		c = out.fit.co;
	else
		c = out.fit.c;
	end		
else
	c = (make_avg_any([s.win.c{:}], navg))*mag_cal;
	
	% Remove outliers
	if(size(c, 1) > 2)
		outliers = arrayfun(@(x)outlier(c(:, x), 0.1), 1:size(c(:, :), 2), ...
								'UniformOutput', false);
							
		outloc = find(~cellfun(@isempty, outliers));
		
		trans = 1:size(c, 1);
		
		if(~isempty(outloc))
			for i = 1:length(outloc)
				ct = trans;
				ct(outliers{outloc}) = [];
				c(outliers{outloc}, outloc) = mean(c(ct, outloc));
			end
		end
		
		
	end
	
	c = squeeze(mean(c, 1));
end

if(c(1) < 0)
	c = -c;
end

out.disp.mag_cal = mag_cal;

% Get the t1 for each measurement.
perm = 1:p.nDims; 
perm(dim) = [];
perm = [dim, perm(:)];

ns = prod(p.maxsteps(1:length(p.maxsteps) ~= dim));

out.fit.co = c;

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

	
out.fit.(name) = zeros(ns, nc);
out.fit.(namea) = zeros(ns, nc);
out.fit.(namee) = zeros(ns, nc);
out.fit.(nameae) = zeros(ns, nc);
out.fit.(nameb) = zeros(ns, nc*2);

out.fit.c = c;
out.fit.cf = zeros(size(c));

out.fit.beta = zeros(ns, length(typical_values));

warning('off'); %#ok;

zm = zeros(size(typical_values));
ub = ones(size(zm))*10;
ub(2:2:end) = 1e10;

for i = 1:ns
	[tau, ~, r, ~, ~, ~, J]  = lsqcurvefit(kern, typical_values, t, c(:, i)', zm, ub, options);

	out.fit.(name)(i, :) = tau(1:2:end);
	out.fit.(namea)(i, :) = tau(2:2:end);
	out.fit.cf(:, i) = kern(tau, t);

	% Standard error (95% confidence interval)
	ci = nlparci(tau, r, 'jacobian', J);
	out.fit.(namee)(i, :) = ci(1:2:end, 2)-tau(1);
	out.fit.(nameae)(i, :) = ci(2, 2)-tau(2);

	out.fit.(nameb)(i, :) = tau;
end


tau = out.fit.(name);
taue = out.fit.(namee);

warning('on'); %#ok