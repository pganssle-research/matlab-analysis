function [out, t1, t, c] = find_t1_lsq(s, mag_cal, navg)
% Adds a t1 to s.fit
%
% Usage:
% [out, t1, t, c] = find_t1_lsq(s, mag_cal, navg);

if(~exist('mag_cal', 'var') || isempty(mag_cal) || mag_cal <= 0)
	if(isfield(s, 'disp') && isfield(s.disp, 'mag_cal'))
		mag_cal = s.disp.mag_cal;
	else
		mag_cal = 836.52;
	end
end

if(~exist('navg', 'var'))
	navg = 6;
end

out = s;

% Find the dimension we're interested in.
p = s.prog;
ins = find(p.vtypes == 1);
dels = [p.vdel{ins}];
[~, i] = max(dels(end, :));
ins = ins(i);

t = p.vdel{ins}'; % Time vector.

c = squeeze(make_avg_any([s.win.ac{:}], navg))'*mag_cal;

if(c(1) < 0)
	c = -c;
end

out.disp.mag_cal = mag_cal;

% Get the t1 for each measurement.
typical_values = [2, c(1)];
perm = 1:p.nDims;
perm(p.vinsdim(ins)) = [];
perm = [p.vinsdim(ins), perm(:)];

ns = prod(p.maxsteps(perm(2:end)));

if(length(perm) > 1)
	c = permute(c, perm); % Put the varying dimension at the front.
end

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
	'Typical', typical_values, 'MaxFunEvals', 2000);
  
% Get the results, save only the exitflag which gives convergence.
out.fit.t = t;
out.fit.c = c;
out.fit.cf = zeros(size(c));
out.fit.t1 = zeros(1, ns);
out.fit.t1A = zeros(1, ns);

out.fit.t1e = zeros(1, ns);
out.fit.t1Ae = zeros(1, ns);

warning('off'); %#ok;

for i = 1:ns
	[tau, ~, r, ~, ~, ~, J]  = lsqcurvefit(@exponential_fit, typical_values, t, c(:, i)', [], [], options);
	
		
	out.fit.t1(i) = t1(1);
	out.fit.t1A(i) = t1(2);
	out.fit.cf(:, i) = exponential_fit(t1, t);

	% Standard error (95% confidence interval)
	ci = nlparci(t1, r, 'jacobian', J);
	out.fit.t1e(i) = ci(1, 2)-t1(1);
	out.fit.t1Ae(i) = ci(2, 2)-t1(2);
end

t1 = out.fit.t1;

warning('on'); %#ok