function fid = SVDSpectrum(in, c, nu)
% This function gives you an FID based on what you give it from linear SVD
% predictions.
%
% nu may not be longer than half the total length of the input
%
% in = data in
% c = number to cut
% nu = number to use
%
% The output will be the same size as the input.
%
% Usage:
% [fid, t] = SVDSpectrum(in, c, nu);

% Pull out the relevant part
data = in(c+1, c+nu);

% Make a Hankel matrix out of it
h1 = hankel(in(c+1, c+nu*2));

h = h1(nu, nu);

% Get the SVD
[u l v] = svd(h);

ld = diag(l);

kdiff = zeros(length(ld)/2, 1);
for j = 1:length(ld)/2
    kdiff(j) = ld(2*j-1)-ld(2*j);
end

% Here we determine if it's leveling out
perc = 0.5; % How much of a change do we need to see before we call it the cutoff
for j = 1:length(kdiff)
    if abs((kdiff(j) - kdiff(j+1))/kdiff(j)) > perc || abs((kdiff(j)-kdiff(j+1))/kdiff(j+1)) > perc
        continue
    else
        break
    end
end

k = j;

% I think we're reshaping the matrices to include the correct rank
ut = u(1:end, 1, k);
wt = w(1:k, 1:k);
vt = v(1:end, 1:k);

% New signal vector I think
newsig = ut*wt*vt';

% Generate something where each element along the antidiagonal is set to be
% the average along the antidiagonals.
sigt = fliplr(newsig); % Make the antidiagonals diagonals

% Know how big this matrix is
s = size(newsig);
s = s(1);

nsigt = zeros(size(sigt));

for n = 1:(2*s-1)
    if n > s
        m = s-n+1;
    else
        m = n-1;
    end
    a = mean(diag(sigt, m));
    
    nsigt(n:s+1:(n+(s+1)*(s-m)) = a;
end








