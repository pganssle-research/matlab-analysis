function s_out = sort_by_title(s_in, descend)
% Given an array of structures with a field 'Title' or 'Name', returns the
% same array of structs, sorted.
%
% If descend evaluates to logical true, the sort will be descending.
% Otherwise it's ascending.
%
% Usage:
% s_out = sort_by_title(s_in, [descend]);

names = {'Title', 'title', 'Name', 'name'};

if(~exist('descend', 'var'))
    descend = false;
end

for j = 1:length(names)
    if(isfield(s_in, names{j}))
        [~, i] = sort(eval(['{s_in.' names{j} '}']));
        break;
    end
end

if(~exist('i', 'var'))
    error('Struct array title could not be found.');
end

if(descend)
    i = fliplr(i);
end

s_out = s_in(i); % Initialize them as the same sort of thing.