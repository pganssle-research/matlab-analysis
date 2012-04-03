function [spec, spec_m, f1, f2] = make_avg_2d_spec(in, t1, t2, cut, lb1, lb2, ds)
% Generates the complex and magnitude spectrum for the data you give it.
%
% in = data, 2 dimensional
% t1 = data vector, direct dimension
% t2 = data vector, indirect dimension
%
% cut = number of points to cut -> default is 0
% lb1, lb2 = line broadening -> default is 0
%
% d2 = downsampling of the direct dimension -> default is 16
%
% Spec = 2D spectrum
% Spec_m = 2D magnitude spectrum
%
%
% Usage: 
%
% [spec, spec_m, f1, f2] = make_avg_2d_spec(in, t1, t2, cut, lb1, lb2, ds);

% Setup the defaults
if nargin < 4
    cut = 0;
end

if nargin < 5
    lb1 = 0;
end

if nargin < 6
    lb2 = 0;
end

if nargin < 7
    ds = 16;
end

si = size(in);
np1 = si(1);
if numel(si) < 3
    np2 = si(2);
else
    in = squeeze(mean(in(:, :, :), 2)); % Get the average
    np2 = si(3);
end

sr1 = np1/t1(end);
sr2 = np2/t2(end);

% First the polynomial subtraction and cut+zero-pack
for n = 1:np2
    in(:, n) = sub_poly(in(:, n), 2, cut); % Polynomial subtraction, 2nd order
    in(:, n) = cutandzeropack(in(:, n), cut); % Cut and zero pack
end

% Now we apodize the data
in = apodize2d(in, t1, t2, lb1, lb2);

% Do the fourier transform on it
newnp1 = 2^(ceil(log2(np1)));
newnp2 = 2^(ceil(log2(np2)));

s = fft2(in, newnp1, newnp2); % The 2D Fourier transform

% If the downsampling is greater than the number of points, don't
% downsample at all.
if newnp1 < ds
    ds = 1;
end

d = newnp1/floor(newnp1/ds); % Make sure it's an even multiple

s_d = s(1:d:end, :); % Downsampled spectrum
si = size(s_d);

newnp1 = si(1); % Need to know the downsampled number of points.

s_d = s_d(1:(newnp1/2), 1:(newnp2/2)); % Just take the first half of each dimension.

% Set the frequency vectors:
f1 = 0:sr1/(newnp1-1):(sr1/2);
f2 = 0:sr2/(newnp2-1):(sr2/2);

% Set the final outputs
spec = s_d;
spec_m = magnitude(spec);
