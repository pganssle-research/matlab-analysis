function [peak_int, ref] = get_peak_int(data_in, ref_start, ref_end)

ref = arrayfun(@(x)(x > ref_start) && (x < ref_end), data_in.x);

peak_int = sum(data_in.data(:, ref), 2);
