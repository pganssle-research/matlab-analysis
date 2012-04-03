function data_out = fix_offsets(data_in, ref_start, ref_end, offset)
% Given the struct data_in of 1D data, this will ensure that all the data
% are properly aligned.
%
% If ref_start and ref_end are 2-vectors, the data must be 2-d, otherwise
% it must be 1-d.

is_1d = logical(~[data_in.is2d]);
is_2d = false;

if(length(ref_start) > 1)
    is_1d = ~is_1d;
    is_2d = true;
end

if(~is_2d)
    data_out = struct('Name', '', 'x',[], 'data', [], 'irefs', []);
else
    data_out = struct('Name', '', 'x', [], 'y', [], 'data', [], 'irefs', []);
end

if(sum(is_1d) == 0)
   data_out = struct();
  return;
end

data_b = data_in(is_1d);
data_out.Name = {data_b.Title};

x = data_b(1).XAxis;    % These need to be the same for this to work
s = num2cell(size(data_b(1).Data));
data = zeros(length(data_b),s{:});

if(is_2d)
    data = permute(data, [1, 3, 2]);
end

for i = 1:length(data_b)
    bd = [data_b(i).Data]';
    data(i, :) = bd(:);
end

refx = arrayfun(@(j)(j > ref_start(1)) && (j < ref_end(1)), x);

if(is_2d)
    y = data_b(1).YAxis;
    refy = arrayfun(@(j)(j > ref_start(2)) && (j < ref_end(2)), y);
    
    % Find the peak in the region    
    [m, ix] = max(data(:, refx, refy), [], 2);
    [~, iy] = max(squeeze(m), [], 2);
    
    ix = ix(iy);
else
    [~, ix] = max(data(:, refx), [], 2);
end

% Shifts to the left, then truncates by the maximum value shifted.
offx = min(ix)-ix;
tnum_x = max(ix)-min(ix);

if(is_2d)
    offy = min(iy)-iy;
    tnum_y = max(iy)-min(iy);
    
    for j = 1:size(data, 1)
        data(j, :, :) = circshift(data(j, :, :), [0, offx(j), offy(j)]);
    end
    
    data(:, :, (end-tnum_y):end) = [];
else
    for j = 1:size(data, 1)
        data(j, :) = circshift(data(j, :), [0, offx(j)]);
    end
end


data(:, (end-tnum_x):end, :) = [];

% Handle the changes in the axis.
x_sub = x(refx);
x((end-tnum_x):end) = [];
x = x - x_sub(min(ix)) + offset(1);
data_out.x = x';
data_out.data = data;

if(is_2d)
    y_sub = y(refy);
    y = y - y_sub(min(iy)) + offset(1);
    
    data_out.y = y';
    data_out.irefs = arrayfun(@(j)sum(abs(x)), data(:, refx, refy));
else
    data_out.irefs = arrayfun(@(j)sum(abs(x)), data(:, refx));
end

data_out.irefs = sum(data_out.irefs(:, :), 2);