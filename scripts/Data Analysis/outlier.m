function [index] = outlier(y, alpha, k)
% Detection and Removal of Outliers in Data Sets
%    ( Rosner's many-outlier test)
%
%       index = outlier( y [, crit, k] )
%
% where  index = indices of outliers in the data
%        y     = data set (should be stationary)
%        crit  = detection criterion (default 2)
%			k     = maximum number of outliers
%
% Originally written by Bob Newell, February 1996
% Modified by Jaco de Groot, May 2006
%
% Bob Newell used a fixed value for lambda. This script calculates the
% critical values for lambda based on the equations in
%
% "Quality control of semi-continuous mobility size-fractionated 
%    particle number concentration data",
% 
% Atmospheric Environment 38 (2004) 3341–3348, Rong Chun Yu,*, 
%     Hee Wen Teh, Peter A. Jaques, Constantinos Sioutas, John R. Froines)
%-----------------------------------------------------
%
% Modified by Paul J. Ganssle, June 2012
% Code is now completely vectorized. I also fixed some errors.

% Vectorize data
y = y(:);
n = length( y ); % Looks like boobs.

% Set up Defaults
if(~exist('alpha', 'var'))
	alpha = 0.05;
end

if(~exist('k', 'var'))
	k = n-2;				% Leave at least 2!
end

% sort deviations from the mean
ybar = mean( y );
[ ys, is ] = sort( abs( y - ybar ));

% Statistics for up to k outliers.
ys(ys == 0) = -NaN; % Just in case there's an actual zero in there.

yy = rot90(triu(repmat(ys, 1, k), k-n), 2)';
yy(yy == 0) = NaN;	% So we can ignore the zeros in the mean and NaN
							% Useful for vectorizing.

yy(yy == -NaN) = 0;
if(isvector(yy))
	R = abs(yy(1) - nanmean(yy, 2))./nanstd(yy, 0, 2);
else
	R = abs(diag(yy)-nanmean(yy, 2))./nanstd(yy, 0, 2);
end


% Statistical test for outliers
n = n-(1:k);
pcrit = 1-(alpha./(2*(n+1)));
t = tinv(pcrit, n-1);
lambda = ((n.*t).*(((t.^2+n-1).*(n+1)).^(-1/2)))';

index = is(n(R>lambda)+1);
