function s = parse_instructions(prog)
% Parse instructions into a meaningful structure. Must pass this something
% with a prog.instrs field.
%
% Inputs:
% prog:	Pulse program type structure with prog.instrs cell array.
%
% Outputs:
% s:		A parsed structure with the cell array broken into structs and
%			full pulse sequences for each point in multidimensional sequences.
%			More machine-readable, less human-readable.
%
% Usage:
% s = parse_instructions(prog);

% This should define things like CONTINUE, STOP, etc.
declare_prog_vars;

instrs = {'CONTINUE', 'STOP', 'LOOP', 'END_LOOP', 'JSR', 'RTS', 'BRANCH', 'LONG_DELAY', 'WAIT'};
u = struct('s', 1, 'ms', 1000, 'us', 1e6, 'ns', 1e9);

p = prog;

s.ni = p.ninst;

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
if(prog.varied && isfield(prog, 'vinslocs'))
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