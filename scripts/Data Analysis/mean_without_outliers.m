function [dmean, dstd, outliers] = mean_without_outliers(data, dim, crit)
% Feed this data and specify along which dimension to find the mean. This
% will remove any outliers and return the mean and standard deviation of
% the data, as well as the number of outliers at each point.
%
% Dimension is 2 if not otherwise specified. If the input is a vector,
% default dimension is the vector one.
%
% data	=	the data, ND data set
% dim		=	dimension along which it varies 
% crit	=	detection criterion (default: 0.1)
%
% [mean, std, outliers] = mean_no_outliers(data, dimension, crit);

if(~exist('dim', 'var'))
	if(iscolumn(data))
		dim = 1;
	else
		dim = 2;
	end
end

if(~exist('crit', 'var'))
	crit = 0.1;
end

s = size(data);
perm = 1:length(s);
perm([1 dim]) = [dim 1];

d = permute(data, perm);
s = size(d);

np = length(d(1, :));

if(size(d, 1) < 3)
	dmean = mean(d, 1);
	dstd = std(d);
	outliers = {};
	
	dmean = reshape(dmean, s(2:end));
else
	outliers = arrayfun(@(x)outlier(d(:, x), crit), 1:np, ...
		'UniformOutput', false);

	if(length(s) == 2)
		s = [s, 1];
	end

	dmean = zeros(s(2:end));
	dstd = zeros(s(2:end));

	t = 1:size(d, 1);
	for i = 1:np
		ct = t;
		ct(outliers{i}) = [];
		dmean(i) = mean(d(ct, i), 1);
		dstd(i) = std(d(ct, i));
	end
end






