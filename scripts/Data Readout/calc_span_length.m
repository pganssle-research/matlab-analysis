function len = calc_span_length(instrs, span)
% Calcualtes the length (in seconds) of a span of instructions. 
%
% Inputs:
% instrs:	An instrs structure as found in mc_struct.prog.ps.instrs
% span:		A 1x2 array of the form [s, e] where s = start of span and e =
%				end of span. Uses a 1-based index. Default is full span.
%
% Outputs:
% len:		Time in seconds.
%
% Usage:
% len = calc_span_length(instrs[, span]);

declare_prog_vars; % Should give values for things like LOOP and LONG_DELAY
is_loop = 0;
len = 0;
spans = [];

if(~exist('span', 'var'))
	span = [1, instrs.ni];
end

if(instrs.instr(span(1)) == LOOP && instrs.instr(span(2)) == END_LOOP && instrs.data(span(2)) == span(1)-1)
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