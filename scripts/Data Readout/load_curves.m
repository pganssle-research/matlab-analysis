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
% .rem_outliers: Boolean, whether or not to remove outliers (default on)
% .off: Offset
%
% [out, t, c] = load_curves(s, options);

if(exist('options', 'var'))
	o = options;
else
	error('Must submit an options structure.');
end

if(~isfield(o, 'type') || isempty(o.type))
	error('Must choose at least one valid type!')
end

if(~isfield(o, 'mag_cal') || isempty(o.mag_cal) || o.mag_cal <= 0)
	if(isfield(s, 'disp') && isfield(s.disp, 'mag_cal'))
		o.mag_cal = s.disp.mag_cal;
	else
		o.mag_cal = 1;
	end
end

if(isfield(o, 'navg') && navg >= 1)
	navg = o.navg;
else
	navg = -1;
end

if(~isfield(o, 'off'))
	o.off = 0.0;
end

if(isfield(o, 'rem_outliers'))
	rem_outliers = logical(o.rem_outliers);
else
	rem_outliers = true;
end

types = {'t1', 't2', 'diff', 'field', 'direct'};
get_type = cellfun(@(t)any(strcmp(t, o.type)), types);
if(~any(get_type))
	error('Must choose a valid type!');
end

is_t1 = logical(get_type(strcmp('t1', types)));
is_t2 = logical(get_type(strcmp('t2', types)));
is_fieldv = logical(get_type(strcmp('field', types)));
is_diff = logical(get_type(strcmp('diff', types)));
is_direct = logical(get_type(strcmp('direct', types)));

if(is_fieldv && (~is_field(o, 'Z_cal') || isempty(o.Z_cal) || o.Z_cal == 0))
	if(isfield(s, 'disp') && isfield(s.disp, 'z_cal') && s.disp.z_cal ~= 0)
		o.Z_cal = s.disp.z_cal;
	else
		o.Z_cal = 1.0;
		warning('Z calibration set to 1.');
	end
end

if(is_diff && (~is_field(o, 'G_cal') || isempty(o.G_cal) || o.G_cal == 0))
	if(isfield(s, 'disp') && isfield(s.disp, 'G_cal') && s.disp.G_cal ~= 0)
		o.G_cal = s.disp.G_cal;
	else
		o.G_cal = 1.0;
		warning('G caibration set to 1.');
	end
end

% Start getting all the curves we need.
has_indirect = any([is_t1, is_t2, is_fieldv, is_diff]);
out = s;

% Get the curves in the direct dimension.
c = s.win.c{:}-o.off/o.mag_cal;
cda = s.win.c{:}*o.mag_cal;
td = s.win.ct{:};

% Need to swap the first two dimensions.
perm = 1:length(size(cda));
perm([2, 1]) = [1, 2];

cda = permute(cda, perm);

if(rem_outliers)
	[cd, stdds] = mean_without_outliers(cda, 1);
else
	cd = mean(cda, 1);
	
	scd = size(cd);
	cd = reshape(cd, scd(2:end));
	stdds = std(cda);
end

stdd = mean(stdds(:));

out.fit.td = td;
out.fit.cd = cd;
out.fit.cda = cda;
out.fit.std = stdd;
out.fit.stds = stdds;

if(has_indirect)
	ci = make_avg_any(cd, navg);
	out.fit.ci = ci;
	
	p = out.prog;
	
	% Get the various indirect mechanisms of variation, if applicable.
	types = p.vtypes(p.vtypes ~= 0);
	
	if(is_t1)
		dim = p.vinsdim(types == 1);
		
		dels = [p.vdel{dim}];
		
		[~, i] = max(dels(end, :));
		
		% If you have identical maxima, hopefully all values are identical, I
		% suppose.
		if(~isscalar(i))
			i = i(1);
		end
		
		dim = dim(i);
		
		t = [p.vdel{dim}]; % Time Vector.
		
		% Rearrange c so that t1 is in the first dimension.
		t1.t = t;
		t1.c = rearrange_ci(ci, dim);
		
		out.fit.t1 = t1;
	end
	
	if(is_t2)
		pins = p.vins(types == 2);
		
		% Find a varying loop
		ins = find(p.ps.instrs.instr(pins+1) == 2, 1);
		
		dim = p.vinsdim(p.vins == pins(ins));
		
		dats = [p.vdata{dim}];
		[~, i] = max(dats(end, :));
		
		if(~isscalar(i))
			i = i(1);
		end
		
		dim = dim(i);
		
		span = find_loop_locs(p.ps, p.vins(ins));
		instrs = p.ps.instrs;
		
		instrs.data(span(1, 1)) = 1;
		
		t2.t = dats*calc_span_length(instrs, span(1, :));
		t2.c = rearrange_ci(ci, dim);
		
		out.fit.t2 = t2;		
	end
end


function c = rearrange_ci(ci, dim)
% Rearranges so that the 'dim' index is the first dimension.
if(isvector(ci))
	if(isrow(ci))
		c = ci';
	else
		c = ci;
	end
else
	perm2 = 1:length(size(ci));
	perm2(dim) = [];
	perm2 = [dim, perm2];
	c = permute(ci, perm2);
end














