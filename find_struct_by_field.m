function [s, loc] = find_struct_by_field(in, name)
% Find a struct from the name of its field.
% Stops at the first one it finds. Breadth-first.

s = [];
loc = [];

flist = fieldnames(in);
num = find(strcmp(name, flist), 1, 'first');

if(~isempty(num))
	s = in.(flist{num});
	loc = [flist{num}];

	return;
end

for i = 1:length(flist)
	b = in.(flist{i});
	if(isstruct(b))
		[s, l] = find_struct_by_field(b, name);
		if(~isempty(s))
			loc = [flist{i} '.' l];
			break;
		end			
	end
end