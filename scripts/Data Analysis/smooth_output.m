function [out, num] = smooth_output(in, alpha)
% Give this an array where the transients are the 2nd dimension and it will
% Go through one-by-one and replace any outliers with the mean value.
%
% in = Data
% stdevs = How many standard deviations is the cutoff (scalar)
% 
% Usage:
% out = smooth_output(in, stdevs);

if(~exist('stdevs', 'var'))
    alpha = 0.1;
end

s = num2cell(size(in));
perm = 1:length(s);
perm([1 2]) = [2 1];

out = permute(in, perm);
s2 = num2cell(size(out));

% Get the mean for each set of transients.
len = length(out(1, :));

locs = zeros(s2{:});
locs(:) = cell2mat(arrayfun(@(x)outlier(out(:, x), alpha), 1:len, 'UniformOutput', false));
locs = logical(locs);

% Recalculate the means based on this
means = zeros(1, s2{2:end});
means(:) = arrayfun(@(x)mean(out(~locs(:, x), x), 1), 1:len);
means = repmat(means, 1, s{1});

% stds(:) = arrayfun(@(x)std(out(:, x)), 1:len); May want to use later.


out(locs) = means(locs);
out = permute(out, perm);
num = sum(locs(:));
