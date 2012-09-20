function ds = bin_peaks(ds, bin_sizes, nDims)
% Provide a cell array of structs with (at least) the fields 'peaks' and
% 'nDims'. This will select all the structs where the nDims field matches
% nDims, then bin the peaks to minimize the number of dimensions.
%
% If nDims is not provided, the mode is used.
% If bin_sizes is a number, an n-square grid of size nxnx(...)xn will be used. 
%
% Usage:
% out = bin_peaks(ds, nDims);

nsets = length(ds);
nDs = zeros(nsets, 1);
peaks = cell(nsets, 1);

i = 1;
for ii = 1:length(ds)
	if isfield(ds{ii}, 'nDims')
		nDs(i) = ds{ii}.nDims;
	else
		continue;
	end
	
	if isfield(ds{ii}, 'peaks')
		peaks{i} = ds{ii}.peaks;
	else
		continue;
	end
	
	i = i+1;
end

nDs(i:end) = [];
peaks(i:end) = [];

if ~exist('nDims', 'var')
	nDims = mode(nDs);
end

% Remove the sets that don't have the right number of dimensions.
peaks(nDs ~= nDims) = [];

% Generate peak and intensity lists now.
intensities = cell(length(peaks), 1);
poses = cell(length(peaks), 1);
for ii = 1:length(peaks)
	pset = peaks{ii};
	pos = pset(:, 1:nDims);
	intensity = pset(:, nDims+1);
	
	poses{ii} = pos;
	intensities{ii} = intensity;
end

% Calculate a grid from the bin sizes.
dinds = cell(nDims, 1);

% Find the minimum and maximum values in the data for each dimension
mind = ones(ii, nDims)*NaN;
maxd = ones(ii, nDims)*NaN;
for ii = 1:length(poses)
	mind(ii, :) = min(poses{ii});
	maxd(ii, :) = max(poses{ii});
end

mind = min(mind);
maxd = max(maxd);

gsize = cell(nDims, 1);
if isscalar(bin_sizes)
	for ii = 1:nDims
		dinds{ii} = linspace(mind(ii), maxd(ii), bin_sizes);
		gsize{ii} = bin_sizes;
	end
else
	for ii = 1:nDims
		dinds{ii} = mind(ii):bin_sizes(ii):maxd(ii);
		gsize{ii} = length(dinds{ii});
	end
end

dgrid = zeros(length(poses), gsize{:});

for ii = 1:length(poses)
	bgrid = zeros(gsize{:});
	pks = poses{ii};
	ints = intensities{ii};
	for jj = 1:length(ints)
		ind = num2cell(ones(1, nDims)*NaN);
		for kk = 1:nDims
			ind{kk} = find(pks(jj, kk) <= dinds{kk}, 1, 'first');
			if pks(jj, kk) < dinds{kk}(ind{kk}) && ind{kk} > 1
				ind{kk} = ind{kk}-1;
			end
		end
		
		bgrid(ind{:}) = bgrid(ind{:}) + ints(jj);
	end
	
	dgrid(ii, :) = bgrid(:);
end

% Remove unused stuff.
ds = dgrid;
ds(:, max(ds) == min(ds)) = [];
return;
	
% Now get the nearest-neighbor values for the positions.
% Don't actually use this - for non-rectangular bins.
if false
	l2dists = cellfun(@(m)L2_distance(m', m'), poses, 'UniformOutput', false);

	for ii = 1:length(l2dists)
		l2dists{ii}(l2dists{ii} == 0) = NaN;
	end

	nns = cellfun(@(m)min(m), l2dists, 'UniformOutput', false);

	nnvals = [];
	for ii = 1:length(nns)
		cnns = [poses{ii}, nns{ii}'];
		nnvals = [nnvals; cnns]; %#ok
	end

	ds = nnvals;
end







