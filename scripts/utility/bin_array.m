function out = bin_array(nums, padding);
% Generates a binary array out of the numbers in nums.
% Optional argument "padding" gives the minimum number of digits, otherwise
% it will be an array of size [length(nums), ceil(log2(max(nums)))+1]

bin_str = dec2bin(nums);
pad = size(bin_str, 2);

if(exist('padding', 'var') && padding > pad)
    pad = padding;
end

out = zeros(length(nums), pad);
s = size(bin_str);
start = pad-s(2);

for y = 1:s(2)
    out(:, y+start) = arrayfun(@(x)str2double(bin_str(x, y)), 1:s(1));
end