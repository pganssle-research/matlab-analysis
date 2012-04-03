function out = 2d_apodize(data, time, lb1, lb2)
% Function that generate apodized data in 2 dimensions
%
%
% data = data you want in
% t1v = time vector for dimension 1 (n x 1)
% t2v = time vector for dimension 2 (n x 1)
% lb1 = line broadening in dimension 1
% lb2 = line broadening in dimension 2
%
% Usage: out = 2d_apodize(data, time, lb1, lb2)

% Generate the exponentials
e_t1 = exp(-t1/lb1);
e_t2 = exp(-t2/lb2);

% Get the outer product to get the 2d version.
app_vec = e_t1*e_t2';

% Apply it to the data
out = app_vec .* data;