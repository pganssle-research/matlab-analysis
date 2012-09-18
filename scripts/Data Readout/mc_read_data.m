function [out, path] = mc_read_data(path, ns)

if(~exist('path', 'var') || isempty(path) || ~ischar(path))
	path = -1;
end

if(~exist('ns', 'var'))
	ns = 1;
end

% Get the raw structure
[s, f] = mc_read_bin(path, 'mc_read_data_hist.mat');

out = [];
if(isempty(f))
	return;
end

% Separately process the groups
MCD_DATAHEADER = '[Data Header]';
MCD_DISPHEADER = '[Display Header]';
MCD_DATAGROUP = '[DataGroup]';

% Data header should come first - That'll be the main portion of the
% structure - so those are top-level values.
s1 = find_struct_by_name(f, MCD_DATAHEADER);
if(isempty(s1))
	return;
end

% Data Structure Names
MCD_FNAME = 'filename';
MCD_ENAME = 'ExperimentName';
MCD_ENUM = 'ExperimentNum';
MCD_DATADESC = 'Description';
MCD_HASH = 'HashCode';
MCD_NCHANS = 'NumChans';
MCD_TSTART = 'TimeStarted';
MCD_TDONE = 'TimeDone';
MCD_CIND = 'CurrentIndex';

out.Filename = [];
out.ExperimentName = [];
out.ExperimentNum = [];
out.hash = [];
out.tstart = [];
out.tdone = [];
out.nc = 0;
out.cind = -1;

sb = find_struct_by_name(s1.data, MCD_FNAME);
if(~isempty(sb))
	fname = deblank(sb.data');
	li = find(fname == '\', 1, 'last');
	if(isempty(li) || li == length(fname))
		out.Filename = fname;
	else
		out.Filename = fname((li+1):end);
	end
end

sb = find_struct_by_name(s1.data, MCD_ENAME);
if(~isempty(sb))
	out.ExperimentName = deblank(sb.data');
end

sb = find_struct_by_name(s1.data, MCD_ENUM);
if(~isempty(sb))
	out.ExperimentNum = sb.data;
end

sb = find_struct_by_name(s1.data, MCD_DATADESC);
if(~isempty(sb))
	out.desc = deblank(sb.data');
end

sb = find_struct_by_name(s1.data, MCD_HASH);
if(~isempty(sb))
	out.hash = sb.data;
end

sb = find_struct_by_name(s1.data, MCD_TSTART);
if(~isempty(sb))
	out.tstart = deblank(sb.data');
end

sb = find_struct_by_name(s1.data, MCD_TDONE);
if(~isempty(sb))
	out.tdone = deblank(sb.data');
end

sb = find_struct_by_name(s1.data, MCD_NCHANS);
if(~isempty(sb))
	out.nc = sb.data;
end

sb = find_struct_by_name(s1.data, MCD_CIND);
if(~isempty(sb))
	out.cind = sb.data;
end

% Display header is next - we can just do a direct dump
out.disp = [];
[~, loc] = find_struct_by_name(f, MCD_DISPHEADER);
if(isfield(s, loc))
	out.disp = eval(['s.' loc]);
end

if(isfield(out, 'disp'))
	% These until these are added to Experimental Parameters
	
	out.disp.z_cal = 0.1379; % Current z calibration
	out.disp.G_cal = 0.0447;  % Current gradient calibration
end

% Read the program.
out.prog = mc_read_prog(s);

if(~isempty(out.prog))
	np = out.prog.np;
	sr = out.prog.sr;
	out.t = linspace(0, np/sr, np);
	
end

% Get the data itself
out.mdata = [];
[~, loc] = find_struct_by_name(f, MCD_DATAGROUP);
if(isfield(s, loc))
	dg = s.(loc);
	if(isstruct(dg))
		fn = fieldnames(dg);
		
		if(length(fn) >= 1)
			if(isfield(out.prog, 'steps'))
				ps = out.prog.steps;
			else
				ps = length(fn);
			end
			
			ci = num2cell(ps);
			data = zeros(length(dg.(fn{1})), ci{:});
			
			% Parse the index - first we need to figure out the number of
			% digits in each number.
			ind = fn{1};
			ns2 = length(ps);
			ndig = (length(ind)-4)/ns2;
						
			for i = 1:length(fn)
				% Now we can parse the actual indexes
				% The +1 is because it's a 1-based index in Matlab, but in C we
				% stored it as a 0-based index.
				ci = num2cell(str2num(reshape(fn{i}(5:end), ndig, ns2)')+1); %#ok;
				data(:, ci{:}) = dg.(fn{i});
			end
					
			out.mdata = data;
			
			if(length(fn) > 1)
				out.adata = mean(out.mdata, 2);
			end
		end
	end
end

if(isfield(out, 'prog') && isfield(out.prog, 'instrs') && out.prog.use_pb)	
	% Find all the loop locations which it could be (starting with the
	% instruction which triggers the scan.
	
	sn = find(out.prog.ps.instrs.scan == 1, 1, 'first');
	tlspans = find_loop_locs(out.prog.ps, sn, 1); 
	
	% For testing for x-y-x-y trains:
	xpulse = bitor(2^6, 2^7);
	ypulse = bitor(2^8, 2^9);
	xy = 0;
		
	if(sn > 0 && ~isempty(tlspans))
		instrs = out.prog.ps.instrs;
		ins = 0;
		r_loop = 0;
		
		for i = 1:size(tlspans, 1)
			% Check for an x-y-x-y train.
			
			if(tlspans(i, 2)-tlspans(i, 1) == 3)
				% Two possibilities here - wait-pulse-wait-pulse or 
				% pulse-wait-pulse-wait
				
				sp = tlspans(i, 1);
				
				if((bitand(xpulse, instrs.flags(sp)) && ...
					bitand(ypulse, instrs.flags(sp+2))) || ...
					(bitand(ypulse, instrs.flags(sp)) && ...
					bitand(xpulse, instrs.flags(sp+2))))
					
					% Pulse-wait-pulse-wait condition
				
					r_loop = 1;
					c_l = [instrs.ts(sp+1), instrs.ts(sp+3)];
					
					cins = [sp+1, sp+3];
					ins = i;	
					xy = 1;
					break;
				elseif((bitand(xpulse, instrs.flags(sp+1)) && ...
					bitand(ypulse, instrs.flags(sp+3))) || ...
					(bitand(ypulse, instrs.flags(sp+1)) && ...
					bitand(xpulse, instrs.flags(sp+3))))
					
					% Wait-pulse-wait-pulse condition
				
					r_loop = 1;
					c_l = [instrs.ts(sp), instrs.ts(sp+2)];
							
					xy = 1;
					cins = [sp, sp+2];
					ins = i;
					break;
				end
			end
			
			for j = tlspans(i, 1):tlspans(i, 2)
				if(instrs.flags(j) == 0 && instrs.ts(j) > 20e-3)
					% Loop located.
					r_loop = 1;
					c_l = instrs.ts(j);
					cins = j;
					ins = i;
					break;
				end
			end
		end
		
		if(r_loop)
			ad = zeros(out.prog.nDims, 1);
			
			if(out.prog.varied && isfield(out.prog, 'vins'))
				vins = out.prog.vins + 1;
				ad1 = arrayfun(@(j)j > sn && j < tlspans(ins,2 ), vins);
				
				for i = 1:out.prog.nDims
					ad(i) = any(ad1(out.prog.vinsdim == i));
				end
				
											
			end
			
			sd = size(out.mdata);
			cc = num2cell(ones(size(sd)));
			
			ad2 = [0; 0; ad];
			ad2 = logical(ad2);
			
			if(sum(ad) == 0)
				ts = 1;
				vinstrs = [];
				ind = [];
			else
				ind = sd(ad2);
								
				ts = prod(ind);
				vinstrs = out.prog.ps.vinstrs;	
			end
			
			out.win.ad = ad2;
			out.win.ind = ind;
			
			ad = logical(ad);
			
			l_l = ones(ts, 1);
			t_t = zeros(ts, 1);
			c_t = zeros(ts, 1);
			e_t = zeros(ts, 1);
					
			if(isfield(out.prog, 'vins') && ~isempty(intersect(cins, out.prog.vins)))
				clb = c_l;
				c_l = zeros(ts, 1);
				c_l(1) = clb;
			else
				c_l = ones(ts, 1)*c_l;
			end
			
			for i = 1:ts
				if(~isempty(vinstrs))
					[cc{ad}] = ind2sub(ind, i);
					instrs = out.prog.ps.vinstrs(cc{:});
				end
				
				l_l(i) = instrs.data(tlspans(ins, 1));
				t_t(i) = calc_span_length(instrs, tlspans(ins, :)); % Get loop length.
				
				if(c_l(i) ~= 0)
					c_l(i, :) = instrs.ts(cins);
				end
				
				c_t(i) = t_t(i)/l_l(i); % Get per-loop length.
				e_t(i) = calc_span_length(instrs, [sn, tlspans(ins,1)-1]);
			end
			
			c_t = t_t./l_l;
					
			c_t = c_t * 1000; % In ms;
		
			% Things to skip.
			ne = 1;
			start = e_t*1000+c_t*ns;
			num_win = floor((out.t(end)*1000 - start)./c_t - ne);

			if(num_win > l_l - ne)
				num_win = l_l - ne;
			end
			
			ecb = 'out.mdata(:, :, ';
		
			cc = num2cell(ones(size(sd)));
			
			out.odata = out.mdata;
			
			for i = 1:ts
				% Generate a command
				if(~isempty(ind))
					[cc{ad2}] = ind2sub(ind, i);
					inds = '';
					for j = 1:length(ad)
						if(~ad(j))
							inds = [inds, ', :'];
						else
							inds = [inds, ', ', num2str(cc{j+2})];
						end
					end

					inds = inds(3:end);

					ec = [ecb, inds, ');'];
					cdata = eval(ec);
				else
					cdata = out.mdata;
				end
				
				asym = 0.9;
				frac = 0.75;
	
				c_l = (c_l*1000)-20; % 20ms of this will be useless.				
				
				if(~xy)
					frac = frac*(c_l./c_t); % Take 75% of remaining fraction.

					[points, out.win.spans{i}, ~, t_c] = get_subset(cdata, ...
						c_t(i), start(i), asym, frac(i), ...
						num_win(i), out.prog.sr);
				else
					% If we're in xy mode, we want to get two subets with
					% different fractions then interpolate them.
					
					cl1 = c_l(1);
					cl2 = c_l(2);
										
					frac1 = frac*(cl1./c_t);
					frac2 = frac*(cl2./c_t);
					
					% This is the tricky part - each one is a fraction of a
					% fraction. We'll assume that the pulse is an insignificant
					% fraction of the total in any given loop.
					%
					% Here's how this works:
					%        --------
					%       |        |       cl1
					%   a1  |        |   a2   |     a3
					% ______|   b1   |________|___________________
					%
					% Normally, asym = a1/(a1+a2), but in this case we need to
					% transform that in such a way that it represents
					% a1/(a1+a2+a3), so as to ignore a3 (the second lobe). We
					% can do this by multiplying by (a1+a2)/(a1+a2+a3).
					%
					% For the second lobe:
					%
					%					              --------
					%						cl1		 |        |
					%		a3 			  |	a1  |        |   a2   
					%___________________| ______|   b1   |________
					%
					% We want to transform it in this diagram such that 
					% asym2 = (a1+a3)/(a1+a2+a3), we can do this by multiplying
					% asym by (a1+a2)/(a1+a2+a3) as before, then adding
					% a3/(a1+a2+a3).
					
					ws1 = frac1*c_t;
					ws2 = frac2*c_t;
					
					cl1 = cl1+20;
					cl2 = cl2+20;
					
					at1 = c_t-ws1;
					at2 = c_t-ws2;
															
					asym1 = asym*(c_t-cl2-ws1)/at1;
					asym2 = asym*(c_t-cl1-ws2)/at2+cl1/at2;
					
					% Do this twice and interpolate later - it may be preferable
					% to change this function to accept multiple asym/frac
					% values to get subsets of subsets.
					[p1, spans1, ~, t_c1] = get_subset(cdata, c_t(i), ...
						start(i), asym1, frac1, ...
						num_win(i), out.prog.sr);
					
					[p2, spans2, ~, t_c2] = get_subset(cdata, c_t(i), ...
						start(i), asym2, frac2, ...
						num_win(i), out.prog.sr);
					
					points = interp_mat(p1, p2, 2);
					out.win.spans{i} = spans1 - spans2;
					
					t_c = interp_mat(t_c1, t_c2);
				end
				points = mean(points, 1);
				s2 = size(points);
				
				if(~isvector(points))
					points = reshape(points, s2(2:end));
				else
					points = points';		% For whatever reason the vectors come out as row vectors.
				end
				
				points = points * (-2*mod(ns+1, 2)+1);
				
				out.win.it{i} = t_c;
	
				out.win.spans{i} = out.win.spans{i} - mean(out.win.spans{i});
				
				t_c = t_c/1000;
			
				if(out.disp.poly_ord >= 0 && out.disp.poly_ord < 99)
					s2 = size(points);
					
					% Get the fits
					warning('off','all');
					pfit = zeros(out.disp.poly_ord+1, size(cdata(:, :), 2));
					ocdata = cdata;
					
					for j = 1:size(cdata(:, :), 2)
						pfit(:, j) = polyfit(t_c, points(:, j), out.disp.poly_ord);
						cdata(:, j) = ocdata(:, j)-polyval(pfit(:, j), out.t');
						points(:, j) = points(:, j)-polyval(pfit(:, j), t_c);
					end
					
					out.win.polyfit{i} = reshape(pfit, [out.disp.poly_ord+1, s2(2:end)]);
					if(~isempty(ind))
						ec = [ecb, inds, ') = cdata;'];
						eval(ec);
					else
						out.mdata = cdata;
					end
					warning('on','all');
				end
				
				if(xy)
					num_win(i) = num_win(i)*2;
				end
				
				c = zeros([num_win(i)-1, s2(2:end)]);
				ct = zeros(num_win(i)-1, 1);
				c(1:2:end, :) = points(1:2:(end-1), :) - points(2:2:end, :);
				c(2:2:end, :) = points(3:2:end, :) - points(2:2:(end-1), :);
				
				ct(1:2:end) = (t_c(1:2:(end-1))+t_c(2:2:end));
				ct(2:2:end) = (t_c(3:2:end)+t_c(2:2:(end-1)));
				
				ct = ct/2;
				
				out.win.p{i} = points;
				out.win.ap{i} = mean(points, 2);
				out.win.c{i} = c;
				out.win.ac{i} = mean(c, 2);
				out.win.ct{i} = ct;
			end
			
			if(~isempty(ind) && length(ind) > 1)
				out.win.c = reshape(out.win.c, ind);	
				out.win.ac = reshape(out.win.ac, ind);
				out.win.ct = reshape(out.win.ct, ind);
				out.win.it = reshape(out.win.it, ind);
				out.win.ap = reshape(out.win.ap, ind);
				out.win.p = reshape(out.win.p, ind);
				out.win.polyfit = reshape(out.win.polyfit, ind);
				out.win.spans = reshape(out.win.spans, ind);
			end
		end
	end
end

if(~isfield(out, 'win'))
	ord = 2;
	
	if(isfield(out, 'disp') && isfield(out.disp, 'polyord'))
		ord = out.disp.polyord;
	end
	
	out = process_data(out, ord);
end
out = add_fft(out);

function [s, loc] = find_struct_by_name(in, name)
% Find a struct from its .name parameter.
s = [];
loc = [];
flist = fieldnames(in);

for i = 1:length(flist)
	b = in.(flist{i});
	
	if(isfield(b, 'name') && strcmp(b.name, name))
		s = b;
		loc = [flist{i}];
		break;
	end
	
	if(isstruct(b.data))
		[s l] = find_struct_by_name(b.data, name);
		if(~isempty(s))
			loc = [flist{i} '.' l];
			break;
		end
	end
end




