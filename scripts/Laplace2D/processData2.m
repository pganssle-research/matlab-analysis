function out = processData2(data, K1, K2, tau1, tau2, varargin)
% Takes data from a struct and two kernel functions, uses them to calculate
% a 2D laplace inversion with Tikhonov regularization.
%
% Paul J. Ganssle, U.C. Berkeley, Pines Lab, June 2012
%
% Released free for use, modification and non-commercial distribution, with
% credit.
%
% Inputs:
% data -> Struct, can be created from make_data_struct()
%   data.x = First dimension (e.g. time in indirect dimension 1) [size nx1]
%   data.y = Second dimension (e.g. time in indirect dimension 2, voltage,
%            etc) [size mx1]
%   data.z = Measurement response [size nxm]
%
% K1, K2 -> Kernel functions. Can be anonymous function or function handles
%           Functions should be of the form K1(x, tau1), K2(y, tau2) OR 
%           {K1(x, tau1, ...), varargin{:}}, etc. 
%
% tau1, tau2 -> The initial vector you'd like to span in your kernel
%               function (this will be cut down severely in the SVD phase).
%
% Options (use 'name', value)
% 'nSVDs' -> Number of SVDs to use in the calculation. Pass either a
%           scalar, which will be used for both kernels, or a 2x1 vector, 
%           where nSVDs(1) is number of SVDs for K1, nSVDs(2) is the number
%           of SVDs for K2. Default is 8. Pass 0 or -1 per dimension for
%           default (e.g. [-1, 10] will give [8, 10]).
%
% 'alpha' -> Choose an alpha for the Tikhonov normalization. Default: 1
%
% Outputs:
% out -> Struct, transformed such that x<->K1, y<->K2
%   out.t1 = First dimension in transformed (K1) space. [size nx1]
%   out.t2 = Second dimension in transformed (K2) space. [size mx1]
%   out.f = Spectrum in the transrmed space [size nxm]
%
% Usage:
% out = processData(data, K1, K2, ... (options));

d = data;

% Process options first.
% Number of SVDs
ns = find(strcmp(varargin, 'nSVDs'), 1, 'first');
svdd = 8;               % Number of SVDs default.
nSVDs = ones(2, 1)*svdd;
if(~isempty(ns))
    ns = varargin{ns+1};
    
    if(~isnumeric(ns))
        error('Number of SVDs must be numeric.');
    else
       if(isscalar(ns))
           nSVDs(:) = ns;
       else
           nSVDs(:) = ns(1:2);
       end
    end
    
    nSVDs(nSVDs <= 0) = svdd;
end

a = find(strcmp(varargin, 'alpha'), 1, 'first');
alpha = 1;
if(~isempty(a))
    a = varargin{a+1};
    if(isnumeric(a) && isscalar(a) && a >= 0)
       alpha = a; 
    else
        warning('Invalid alpha - using default value of %d.', alpha); %#ok
    end
end

% Generate the matrices for the kernel functions.
if(iscell(K1))
    k1 = K1{1};
    k1 = k1(d.x, tau1, K1{2:end});
else
    k1 = K1(d.x, tau1);
end

if(iscell(K2))
    k2 = K2{1};
    k2 = k2(d.y, tau2, K2{2:end});
else
    k2 = K2(d.y, tau2);
end


% Get the SVDs now:
[U1, S1, D1] = svds(k1, nSVDs(1));
[U2, S2, D2] = svds(k2, nSVDs(2));

% Do the minimization.
Mc = U1'*d.z*U2;
F0 = ones(nSVDs(1), nSVDs(2));

alpha = 2;

o = optimset('Algorithm', 'sqp', 'GradObj', 'on', 'Hessian',  'on', 'MaxIter', 5000, ...
	'TolX', 1e-9, 'LargeScale', 'on', 'TolFun', ...
    1e-16, 'MaxFunEvals', 5000);

K1r = S1*D1';
K2r = S2*D2';

K = K1r*K2r';
K = K(:)';

t = cputime;

for i = 1:9
	alpha = 10^(i-3);

	F = fmincon(@(F)fun(F, Mc, K1r, K2r, alpha), F0, -K, 0, [], [], [], [],[], o);

	F0 = F;
	
	subplot(3, 3, i);
	
	image(tau1, tau2, K1r'*F*K2r, 'CDataMapping', 'scaled');
end
fprintf('Elapsed time %d minutes.\n', (cputime-t)/60);


out.t1 = tau1;
out.t2 = tau2;
out.f = K1r'*F*K2r;
end

function [f, g, H] = fun(F, Mc, K1r, K2r, alpha)
	f = norm(Mc-F, 'fro').^2 + alpha*norm(K1r'*F*K2r, 'fro').^2;
	
	if(nargout > 1)
		g = 2*(Mc-1)*norm(Mc-F, 'fro') + 2*alpha*K1r*K2r'*norm(F, 'fro');
	end
	
	if(nargout > 2)
		H = 2*(Mc);
	end
end