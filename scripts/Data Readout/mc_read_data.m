function [out, path] = mc_read_data(path)

if(~exist('path', 'var'))
	path = -1;
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
MCD_PROGHEADER = '[PulseProgram]';

MCD_NDPC = 'NDPC';
MCD_ANALOGOUT = 'AnalogOutput';
MCD_PULSEPROPS= 'Properties';
MCD_INSTRUCTIONS = 'Instructions';

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

% Pulse program names
MCD_STEPS = 'steps';
MCD_MAXSTEPS = 'maxsteps';
MCD_DATAEXPRS = 'dataexprs';
MCD_DELAYEXPRS = 'delayexprs';
MCD_VINS = 'vins';
MCD_VINSDIM = 'vinsdim';
MCD_VINSMODE = 'vinsmode';
MCD_VINSLOCS = 'vinslocs';

MCD_AOVARIED = 'aovaried';
MCD_AODIM = 'aodim';
MCD_AOVALS = 'aovals';

CONTINUE = 0;
STOP = 1;
LOOP = 2;
END_LOOP = 3;
BRANCH = 6;
LONG_DELAY = 7;
WAIT = 8;

out.FileName = [];
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
		out.FileName = fname;
	else
		out.FileName = fname((li+1):end);
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

% Pulse program
% TODO: Generalize
out.prog = [];
[~, loc] = find_struct_by_name(f, MCD_PROGHEADER);
if(isfield(s, loc))
	out.prog = s.(loc).Properties;
	
	if(isfield(s.(loc), MCD_NDPC))
		s1 = s.(loc).(MCD_NDPC);
		if(isfield(s1, MCD_MAXSTEPS))
			out.prog.maxsteps = s1.(MCD_MAXSTEPS);
		end
		
		if(isfield(s1, MCD_STEPS))
			out.prog.steps = s1.(MCD_STEPS);
		end
		
		if(isfield(s1, MCD_VINS))
			out.prog.vins = s1.(MCD_VINS);
		end
		
		if(isfield(s1, MCD_VINSDIM))
			out.prog.vinsdim = s1.(MCD_VINSDIM);
		end
		
		if(isfield(s1, MCD_VINSMODE))
			out.prog.vinsmode = s1.(MCD_VINSMODE);
		end
		
		if(isfield(s1, MCD_VINSLOCS))
			out.prog.vinslocs = s1.(MCD_VINSLOCS);
		end
		
		if(isfield(s1, MCD_DELAYEXPRS))
			out.prog.delayexprs = s1.(MCD_DELAYEXPRS);
		end
		
		if(isfield(s1, MCD_DATAEXPRS))
			out.prog.dataexprs = s1.(MCD_DATAEXPRS);
		end
		
		
	end
	
	if(isfield(s.(loc), MCD_ANALOGOUT))
		s1 = s.(loc).(MCD_ANALOGOUT);
		if(isfield(s1, MCD_AOVALS))
			out.prog.aovals = s1.(MCD_AOVALS);
		end
		
		if(isfield(s1, MCD_AOVARIED))
			out.prog.aovaried = s1.(MCD_AOVARIED);
			
			if(out.prog.aovaried && isfield(s1, MCD_AODIM))
				out.prog.aodim = s1.(MCD_AODIM);
			end
		end
	end
	
	if(isfield(s.(loc), MCD_INSTRUCTIONS))
		s1 = uint8((s.(loc).(MCD_INSTRUCTIONS))');
		nfields = typecast(s1(1:4), 'int32');
		out.prog.instrs = cell(out.prog.nUniqueInstrs+1, nfields);
		sizes = zeros(nfields, 1);
		types = cell(nfields, 1);
		
		j=5;
		for i = 1:nfields
			l = typecast(s1(j:j+3), 'int32'); % Get the length of the field name
			if(l > 10000)
				error('Memory overload.');
			end
			
			j = j+4;
			
			% Flags
			out.prog.instrs{1, i} = deblank(char(s1(j:j+l-2)));
			j = j+l;
			
			if(strncmp(out.prog.instrs{1, i}, 'instr_data', length('instr_data')))
				out.prog.instrs{1, i} = 'data';
			end
			
			if(strncmp(out.prog.instrs{1, i}, 'trigger_scan', length('trigger_scan')))
				out.prog.instrs{1, i} = 'scan';
			end
			
			if(strncmp(out.prog.instrs{1, i}, 'instr_time', length('instr_time')))
				out.prog.instrs{1, i} = 'time';
			end
			
			if(strncmp(out.prog.instrs{1, i}, 'time_units', length('time_units')))
				out.prog.instrs{1, i} = 'units';
			end
			
			type = typecast(s1(j), 'uint8');
			sizes(i) = fs_size(type);
			types{i} = fs_type(type);
			
			j = j+1;
		end
		
		
		units = {'ns', 'us', 'ms', 's'};
		for i=1:out.prog.nUniqueInstrs
			for k=1:nfields
				out.prog.instrs{i+1, k} = typecast(s1(j:(j+sizes(k)-1)), types{k});
					
				j = j+sizes(k);
			end
			
			out.prog.instrs{i+1, 5} = out.prog.instrs{i+1, 5}*10^(-double(out.prog.instrs{i+1, 6})*3);
			out.prog.instrs{i+1, 6} = units{out.prog.instrs{i+1, 6}+1};
		end
	end
end


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
			
			for i = 1:length(fn)
				[ci{:}] = ind2sub(ps', i);
				data(:, ci{:}) = dg.(fn{i});
			end
			
			out.mdata = data;
			
			if(length(fn) > 1)
				out.adata = mean(out.mdata, 2);
			end
		end
	end
end

if(isfield(out, 'prog') && isfield(out.prog, 'instrs'))
	out.prog.ps = parse_instructions(out.prog);
	
	spans = find_loop_locs(out.prog.ps);
	
	ni = out.prog.ps.ni;
	
	sn = find(out.prog.ps.instrs.scan == 1, 1, 'first');
	
	tlspans = spans;
	a = zeros(size(spans, 1), 1);
	for i = 1:size(spans, 1)
		a(i) = ~isempty(find(arrayfun(@(x, y)spans(i, 1) > x && spans(i, 1) < y, spans(:, 1), spans(:, 2)), 1));
	end
	
	tlspans(logical(a), :) = [];
	
	if(sn > 0 && ~isempty(tlspans))
		instrs = out.prog.ps.instrs;
		ins = 0;
		
		for i = 1:size(tlspans, 1)
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
			
			if(out.prog.varied)
				for j = 1:length(out.prog.vins)
					if(~isempty(find(arrayfun(@(x, y) out.prog.vins(j) >= x && out.prog.vins(j) <= y, spans(:, 1), spans(:, 2)), 1)))
						% Need to recapitulate along this dimenson.
						ad(out.prog.vinsdim(j)) = 1;
					end
				end
			end
			
			sd = size(out.mdata);
			cc = num2cell(ones(size(sd)));
			
			ad2 = [0; 0; ad];
			ad2 = logical(ad2);
			
			if(sum(ad) == 0)
				ts = 1;
				vinstrs = [];
			else
				ind = sd(ad2);
								
				ts = prod(ind);
				vinstrs = out.prog.ps.vinstrs;	
			end
			
			ad = logical(ad);
			
			l_l = zeros(ts, 1);
			t_t = zeros(ts, 1);
			c_t = zeros(ts, 1);
			e_t = zeros(ts, 1);
			
			if(~isempty(find(cins == out.prog.vins, 1)))
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
					c_l(i) = instrs.ts(cins);
				end
				
				c_t(i) = t_t(i)/l_l(i); % Get per-loop length.
				e_t(i) = calc_span_length(instrs, [sn, tlspans(ins,1)-1]);
			
			end
			
			c_t = t_t./l_l;
					
			c_t = c_t * 1000; % In ms;
			c_l = (c_l*1000)-20; % 20ms of this will be useless.				
						
			frac = 0.75*(c_l./c_t); % Take 75% of remaining fraction.
			asym = 0.85;
			
			% Things to skip.
			ns = 2;
			ne = 5;
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
				
				[points, out.win.spans{i}, ~, t_c] = get_subset(cdata, c_t(i), start(i), asym, frac(i), num_win(i), out.prog.sr);
				points = mean(points, 1);
				s2 = size(points);
				points = reshape(points, s2(2:end));
				out.win.it{i} = t_c;
	
				out.win.spans{i} = out.win.spans{i} - mean(out.win.spans{i});
				
				t_c = t_c/1000;
			
				if(out.disp.polyord >= 0 && out.disp.polyord < 99)
					s2 = size(points);
					% Get the fits
					warning('off','all');
					pfit = zeros(out.disp.polyord+1, size(cdata(:, :), 2));
					ocdata = cdata;
					
					for j = 1:size(cdata(:, :), 2)
						pfit(:, j) = polyfit(t_c, points(:, j), out.disp.polyord);
						cdata(:, j) = ocdata(:, j)-polyval(pfit(:, j), out.t');
					end
					
					out.win.polyfit{i} = reshape(pfit, [out.disp.polyord+1, s2(2:end)]);
					ec = [ecb, inds, ') = cdata;'];
					eval(ec);
					
					warning('on','all');
				end
				
				c = zeros([num_win(i)-1, s2(2:end)]);
				ct = zeros(num_win(i)-1, 1);
				c(1:2:end, :) = points(1:2:(end-1), :) - points(2:2:end, :);
				c(2:2:end, :) = points(3:2:end, :) - points(2:2:(end-1), :);
				
				ct(1:2:end) = (t_c(1:2:(end-1))+t_c(2:2:end));
				ct(2:2:end) = (t_c(3:2:end)+t_c(2:2:(end-1)));
				
				ct = ct/2;
				
				out.win.c{i} = c;
				out.win.ac{i} = mean(c, 2);
				out.win.ct{i} = ct;
			end
		end
	end
end

%out = add_fft(out, 2);

function o = fs_type(type)
% File types
FS_CHAR = 1;
FS_UCHAR = 2;
FS_INT = 3;
FS_UINT = 4;
FS_FLOAT = 5;
FS_DOUBLE = 6;
FS_INT64 = 7;
FS_UINT64 = 8;

if(type < 1 || type > 8 || type == FS_CHAR)
	o = 'char';
elseif(type == FS_UCHAR)
	o = 'int8';
elseif(type == FS_INT)
	o = 'int32';
elseif(type == FS_UINT)
	o = 'uint32';
elseif(type == FS_FLOAT)
	o = 'float';
elseif(type == FS_DOUBLE)
	o = 'double';
elseif(type == FS_INT64)
	o = 'int64';
elseif(type == FS_UINT64)
	o = 'uint64';
end

function o = fs_size(type)
% File types
FS_CHAR = 1;
FS_UCHAR = 2;
FS_INT = 3;
FS_UINT = 4;
FS_FLOAT = 5;
FS_DOUBLE = 6;
FS_INT64 = 7;
FS_UINT64 = 8;

if(type < 1 || type > 8)
	o = 1;
elseif(type == FS_CHAR || type == FS_UCHAR)
	o = 1;
elseif(type == FS_INT || type == FS_UINT || type == FS_FLOAT)
	o = 4;
elseif(type == FS_DOUBLE || type == FS_INT64 || type == FS_UINT64)
	o = 8;
end


function [s, loc] = find_struct_by_name(in, name)
% Find a struct from its .name parameter.
s = [];
loc = [];
flist = fieldnames(in);

for i = 1:length(flist)
	b = eval(['in.' flist{i} ';']);
	
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

function s = parse_instructions(prog)
% Parse instructions into a meaningful structure.

CONTINUE = 0;
STOP = 1;
LOOP = 2;
END_LOOP = 3;
JSR = 4;
RTS = 5;
BRANCH = 6;
LONG_DELAY = 7;
WAIT = 8;

instrs = {'CONTINUE', 'STOP', 'LOOP', 'END_LOOP', 'JSR', 'RTS', 'BRANCH', 'LONG_DELAY', 'WAIT'};
u = struct('s', 1, 'ms', 1000, 'us', 1e6, 'ns', 1e9);

p = prog;

s.ni = p.ninst;
s.nUI = prog.nUniqueInstrs;
nUI = s.nUI;

ib = zeros(s.ni, 1);
cb = {cell(s.ni, 1)};
cprog = struct('ni', s.ni, 'tot_time', 0, 'flags', ib, 'instr', ib, 'data', ib, 'time', ib, 'units', cb, 'ts', ib, 'un', ib, 'instr_txt', cb);

% Parse the first version of the program
for i = 2:(s.ni+1)
	units = p.instrs{i, 6};
	un = u.(units);
	time = p.instrs{i, 5};
	ts = time/un;
	
	instr = p.instrs{i, 2};
	data = p.instrs{i, 3};
	
	j = i-1;
	
	cprog.flags(j) = p.instrs{i, 1};
	cprog.scan(j) = p.instrs{i, 4};
	cprog.instr(j) = instr;
	cprog.instr_txt{j} = instrs{instr+1};
	cprog.data(j) = data;
	cprog.time(j) = time;
	cprog.units{j} = {units};
	cprog.un(j) = un;
	cprog.ts(j) = ts;
end

% spans = find_loop_locs(cprog);
s.instrs = cprog;

% Generate a set of pulse program instructions for each step in the
% multi-dimensional space.
if(prog.varied)
	msteps = num2cell(p.maxsteps);
	s.msteps = msteps;
	p.vInstrs = repmat(cprog, msteps{:});
	vil = reshape(p.vinslocs, p.maxnsteps, p.nVaried);
	nv = p.nVaried;
	nis = p.maxnsteps;
	
	cs = msteps;
	for i = 1:nis
		[cs{:}] = ind2sub(p.maxsteps, i);
		for j = 1:nv
			k = vil(i, j)+2; % +1 for non-zero index, +1 for header.
			l = p.vins(j)+1;
			p.vInstrs(cs{:}).flags(l) = p.instrs{k, 1};
			p.vInstrs(cs{:}).instr(l) = p.instrs{k, 2};
			p.vInstrs(cs{:}).data(l) = p.instrs{k, 3};
			p.vInstrs(cs{:}).scan(l) = p.instrs{k, 4};
			p.vInstrs(cs{:}).time(l) = p.instrs{k, 5};
			p.vInstrs(cs{:}).units{l} = p.instrs{k, 6};
			units = p.instrs{k, 6};
			time = p.instrs{k, 5};
			un = u.(units);
			
			p.vInstrs(cs{:}).un(l) = un;
			p.vInstrs(cs{:}).ts(l) = time/un;
		end
		
		s.vinstrs = p.vInstrs;
	end
	
end


function o = data_lims(instr, data)
LOOP = 2;
END_LOOP = 3;
JSR = 4;
RTS = 5;
BRANCH = 6;
LONG_DELAY = 7;

o = 1;
if((instr == LOOP || instr == LONG_DELAY) && data < 2)
	o = 0;
end

o = logical(o);

function len = calc_span_length(instrs, span)
LOOP = 2;
END_LOOP = 3;
LONG_DELAY = 7;
is_loop = 0;

len = 0;
spans = [];

if(instrs.instr(span(1)) == 2 && instrs.instr(span(2)) == 3 && instrs.data(span(2)) == span(1)-1)
	spans = find_loop_locs(instrs, [span(1)+1, span(2)-1], 1);
	l_dat = instrs.data(span(1));
	is_loop = 1;
end

for i = span(1):span(2)
	if(~isempty(spans) && ~isempty(find(arrayfun(@(x, y)i >= x && i <= y, spans(:, 1), spans(:, 2)), 1)))
		continue;
	end

	if(instrs.instr(i) == LONG_DELAY)
		len = len + instrs.ts(i)*instrs.data(i);
	else
		len = len + instrs.ts(i);
	end
end

for i = 1:size(spans, 1)
	len = len + calc_span_length(instrs, spans(i, :));
end

if(is_loop)
	len = len * l_dat;
end

function spans = find_loop_locs(instrs, span, top_level)
if(~exist('top_level', 'var'))
	top_level = 0;
end

if(~exist('span', 'var'))
	span = [1, instrs.ni];
end

LOOP = 2;
END_LOOP = 3;

spans = [];

i = span(1);

ni = span(2);

if(isfield(instrs, 'instrs'))
	instrs = instrs.instrs;
end

while i <= span(2)
	if(instrs.instr(i) == LOOP)
		for j = (i+1):ni
			if(instrs.instr(j) == END_LOOP && instrs.data(j) == i-1)
				% Found it
				spans = [spans; i, j]; %#ok
				
				if(top_level)
					i = j;
				end
			end
		end
	end
	
	i = i+1;
end




