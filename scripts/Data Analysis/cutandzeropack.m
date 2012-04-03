function out = cutandzeropack(in, num)
% Feed this an nd matrix and the number of points you want to take out from
% the beginning. It will return a 2d matrix of the same size with zeros at
% the end.
%
% Usage:
% out = cutandzeropack2d(in, num);

% Get the size for later use
s = num2cell(size(in));

% If they give you the wrong input, just don't truncate at all
if(num < 1)
    num = 0;
end

% First the cut
out = cell2mat(arrayfun(@(x)in((num+1):end, x), 1:length(in(1, :)), 'UniformOutput', false));

% Now the zero pack;
out = reshape(out, s{1}-num, s{2:end});
out = cat(1,  out, zeros(num, s{2:end}));