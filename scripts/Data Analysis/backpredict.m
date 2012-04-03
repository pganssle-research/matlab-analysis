function out = backpredict(in, cut, num_used)
% This cuts out the first 'cut' points, then uses the next 'num_used' points
% in a singular value decomposition to predict what those first points
% should have been.
%
%
% in = data in
% cut = number of points to cut
% num_used = number of points to use for the prediction
%
% Usage:
% out = backpredict(in, cut, num_used);

nu = num_used; % More convenient
c = cut; % More convenient

% Get the data
data = in(c+1:c+nu);

% Now we iteratively generate the points.
points = zeros(c, 1); % Preallocate the points
new_data = data;

for n = 1:c
   % Generate the Hankel matrix
   h = hankel(new_data);

   % Do the SVD
   [u, l, v] = svd(h);

   % Create the q operator:
   q_op = v\l*u';

   q = q_op*new_data;
   q = q/norm(q);
   
   points(c-n+1) = sum(q.*new_data);
   new_data(2:end) = new_data(1:end-1);
   new_data(1) = points(c-n+1);
end

out = in;
out(1:c) = points;