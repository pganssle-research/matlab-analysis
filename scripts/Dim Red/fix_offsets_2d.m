function [data_out, peak_int] = fix_offsets_2d(data_in, ref_start, ref_end, offset)
% Given the struct data_in of 1D data, this will ensure that all the data
% are properly aligned.

is_1d = logical(~[data_in.is2d]);
x = data_in(1).XAxis;
data = [data_in(is_1d).Data]';

ref = arrayfun(@(y)(y > ref_start) && (y < ref_end), x);
[~, i] = max(data(:, ref), [], 2);

% We'll move this stuff to the left and truncate that way.
i = i';
off = min(i) - i;
t_num = max(i)-min(i);

% Shifts to the left, then truncates by the maximum value shifted.
for j = 1:length(i)
    data(j, :) = circshift(data(j, :), [0, off(j)]);
end
data(:, (end-t_num):end) = [];

x_sub = x(ref);
peak_int = sum(data(:, ref), 2);

x = x - x_sub(min(i)) + offset;
x((end-t_num):end) = [];



data_out.Name = {data_in.Title};
data_out.x = x';
data_out.data = data;

