function [out, t2, t, c] = find_t2_lsq(s, mag_cal, navg)
% Adds a t1 to s.fit
%
% Usage:
% [out, t1, t, c] = find_t1_lsq(s, mag_cal, navg);

if(~exist('mag_cal', 'var') || isempty(mag_cal) || mag_cal <= 0)
	if(isfield(s, 'disp') && isfield(s.disp, 'mag_cal'))
		mag_cal = s.disp.mag_cal;
	else
		mag_cal = 893;
	end
end

if(~exist('navg', 'var'))
	navg = 6;
end

out = s;

% Find the dimension we're interested in.
p = s.prog;
ins = find(p.ps.instrs.instr(p.vins(p.vtypes == 2)+1) == 2);
dats = [p.vdata{ins}];
[~, i] = max(dats(end, :));
ins = ins(i);

t = p.vdata{ins}'; % Time vector.

% Now convert it into time
span = find_loop_locs(p.ps, ins);
instrs = p.ps.instrs;
instrs.data(span(1, 1)) = 1;

tau = calc_span_length(instrs, span(1, :));

t = t*tau;


c = squeeze(make_avg_any([s.win.ac{:}], navg))'*mag_cal;

if(c(1) < 0)
	c = -c;
end

out.disp.mag_cal = mag_cal;

% Get the t1 for each measurement.
typical_values = [2, c(1)];
options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'Typical', typical_values, 'MaxFunEvals', 2000);

  
% Get the results, save only the exitflag which gives convergence.
warning('off'); %#ok;
t2 = lsqcurvefit(@exponential_fit, typical_values, t, c, [], [], options);
warning('on'); %#ok

out.fit.t = t;
out.fit.c = c;
out.fit.cf = exponential_fit(t2, t);
out.fit.t2 = t2(1);
out.fit.t2A = t2(2);