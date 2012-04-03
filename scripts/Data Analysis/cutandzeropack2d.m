function out = cutandzeropack2d(in, num)
% Feed this a 2d matrix and the number of points you want to take out from
% the beginning. It will return a 2d matrix of the same size with zeros at
% the end.
%
% Usage:
% out = cutandzeropack2d(in, num);

% Get the size for later use
s = size(in);
np = s(1);

o1 = zeros(s);

% If they give you the wrong input, just don't truncate at all
if(num < 1)
    num = 1;
end

% Do the main function stuff
o1(1:(np-num+1), :) = in(num:end, :);
o1((np-num+1):end, :) = 0;


% Let the output work
out = o1;