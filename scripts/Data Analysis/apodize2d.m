function out = apodize2d(data, t1, t2, lb1, lb2)
% Function that generate apodized data in 2 dimensions
%
%
% data = data you want in
% t1 = time vector for dimension 1 (n x 1)
% t2 = time vector for dimension 2 (n x 1)
% lb1 = line broadening in dimension 1
% lb2 = line broadening in dimension 2
%
% Usage: out = apodize2d(data, time, lb1, lb2)

% Generate the exponentials
e_t1 = exp(-t1/lb1);
e_t2 = exp(-t2/lb2);

% Get the outer product to get the 2d version.
s1 = size(e_t1);
s2 = size(e_t2);

% Make sure that you actually are getting the outer product
% e_t1 needs to be n x 1
if s1(1) == 1
    e_t1 = transpose(e_t1);
end

% e_t2 needs to be 1 x n
if s2(2) == 1
    e_t2 = transpose(e_t2);
end

app_vec = e_t1 * e_t2;

% Apply it to the data
out = app_vec .* squeeze(data);