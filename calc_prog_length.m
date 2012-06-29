function l = calc_prog_length(path_or_prog, verbose, nt)
% Pass this either a program structure, a path, or nothing, and it will
% calculate the length in hours, minutes and seconds that the program will
% take to execute.
%
% Returns the length in seconds. If verbose = 1, prints the length in
% hours, minutes and seconds.
%
% Pass the optional argument nt to see the length with a different number
% of transients. Otherwise the value in the program is used.
%
% Usage;
% calc_prog_length([path_or_prog, verbose, nt]);

histpath = 'hist_calc_prog_len.mat';

if(~exist('path_or_prog', 'var') || isempty(path_or_prog) || ...
		~(isstruct(path_or_prog) || ischar(path_or_prog)))
	path = get_path(histpath, {'.mcd', '.pp'});
elseif(isstruct(path_or_prog))
	prog = path_or_prog;
elseif(ischar(path_or_prog))
	path = path_or_prog;
end

if(~exist('verbose', 'var') || isempty(verbose))
	verbose = 1;
end

if(exist('path', 'var') && ischar(path))
	prog = mc_read_prog(path);
end

if(isfield(prog, 'prog'))
	prog = prog.prog;
end

if(~isfield(prog, 'ps'))
	prog.ps = parse_instructions(prog);
end

p = prog;

if(~exist('nt', 'var'))
	if(isfield(p, 'nt'))
		nt = p.nt;
	else
		nt = 1;
	end
end

l = 0;
if(isfield(p.ps, 'vinstrs'))
	for i = 1:length(p.ps.vinstrs(:))
		l = l+calc_span_length(p.ps.vinstrs(i));
	end
else
	l = calc_span_length(p.ps.instrs);
end

l = l*nt;

if(verbose)
	t = l;
	day = floor(t/86400); t = mod(t, 86400);
	hr = floor(t/3600); t = mod(t, 3600);
	mi = floor(t/60); t = mod(t, 60);
	se = floor(t);

	form = ['Total length of program is:'];
	vars = {};
	
	if(day > 0)
		form = [form, ' %g days '];
		vars = [vars, {day}];
	end
	
	if(day > 0 || hr  > 0)
		form = [form, ' %02g hours '];
		vars = [vars, {hr}];
	end
	
	form = [form, '%02g minutes %02g seconds.\n'];
	vars = [vars, {mi, se}];
	
	fprintf(form, vars{:});
end
