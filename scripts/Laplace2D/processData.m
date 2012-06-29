function out = processData(data, K1, K2, tau1, tau2, varargin)
% This function takes a data structure (x, y, z) and calculates its 2D
% Laplace inversion with Tikhonov inversion. This is an ill-posed problem,
% and boils down to executing the minimization:
%
%	argmin[F] (||Z - K1FK2'||^2 + a*||F||)
%
% As this could be a VERY time-consuming process with large data-sets, the
% data are first compressed to include only the singular values of the
% kernel functions with the largest probability density.
%
% Paul J. Ganssle, U.C. Berkeley, Pines Lab, June 2012
% This function largely implements the algorithm detailed in the paper:
%
% Venkataramanan, L., Song, Y.-Q., Hurlimann, M, D.,
%	"Solving Fredholm Integrals of the First Kind with Tensor Product
%		Structure in 2 and 2.5 Dimensions",
%
%		IEEE Transactions on Signal Processing, Vol. 50, No. 5 May 2002
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
% 'nSVDs' -> Maximum number of SVDs to use in this calculation. Pass either
%				a scalar, which will be used for both kernels, or a 2x1 vector,
%           where nSVDs(1) is number of SVDs for K1, nSVDs(2) is the number
%           of SVDs for K2. Default is [length(data.x), length(data.y)].
%				Pass 0 or -1 per dimension for default.
%
% 'nAlphas' -> Number of alphas to try out Default 16 (best to choose a
%					perfect square).
%
% 'alpha'	-> Initial alpha choice (default 1)
%
% 'alpha_conv' -> Condition for alpha convergence (as a percentage from
%						optimum
%
% 'xscale' -> Pass either 'log' or 'linear', default: linear
%
% 'yscale' -> Pass either 'log' or 'linear' default: linear
%
% 'xlabel' -> Label for the x axis (no default)
%
% 'ylabel' -> Label for the y axis (no default)
%
% 'verbose' -> Boolean - if you want printouts of progress. Default: true
%
% 'dataperc' -> Lower limit on the value of diag(S1)*diag(S2)'. Default:
%					1e-4
%
% 'guess'	-> An initial guess for the spectrum - should be size
%					[length(tau1), length(tau2)].
%
% 'optims'   -> Pass an optimset that will be (indirectly) passed to the
%					 fmincon. Default values set by this function are:
%
%					'GradObj' = 'on'
%					'Hessian' = 'on'
%					'TolX' = 1e-9
%					'TolFun' = 1e-12
%					'Display', 'off'
%
%
% Outputs:
% out -> Struct, transformed such that x<->K1, y<->K2
%   out.t1 = First dimension in transformed (K1) space. [size nx1]
%   out.t2 = Second dimension in transformed (K2) space. [size mx1]
%   out.f = Spectrum in the transrmed space [size nxm]
%   out.F = The spectrum in the compressed domain.
%   out.K = The kernel functions (combined) in the compressed domain.
%
% Usage:
% out = processData(data, K1, K2, ... (options));

d = data;

% Process options first.
% Max number of SVDs - default is max num possible.
ns = find(strcmp(varargin, 'nSVDs'));
svdd = [length(d.x), length(d.y)];
nSVDs = svdd;
if(~isempty(ns))
	NS = varargin{ns(1)+1};
	
	if(~isnumeric(NS))
		error('Number of SVDs must be numeric.');
	else
		if(isscalar(NS))
			nSVDs(:) = NS;
		else
			nSVDs(:) = NS(1:2);
		end
	end
	
	nSVDs(nSVDs <= 0) = svdd(nSVDs <= 0);

	varargin([ns, ns+1]) = [];
end

% Number of alphas
a = find(strcmp(varargin, 'nAlphas'));
nAlphas = 16;
if(~isempty(a))
	na = varargin{a(1)+1};
	if(isnumeric(na) && isscalar(na) && na > 0)
		nAlphas = na;
	else
		warning('Invalid nAlphas - using default value of %d.', nAlphas); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% Initial Alpha
a = find(strcmp(varargin, 'alpha'), 1, 'first');
alpha = 1;
if(~isempty(a))
	al = varargin{a(1)+1};
	if(isscalar(al) && isnumeric(al) && al > 0)
		alpha = al;
	else
		warning('Invalid alpha guess, using default value of %d.', alpha); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% Alpha convergence level
a = find(strcmp(varargin, 'alpha_conv'));
alpha_conv = 0.1; % 10% is a pretty good start.
if(~isempty(a))
	na = varargin{a(1)+1};
	if(isnumeric(na) && isscalar(na) && na > 0)
		alpha_conv = na;
	else
		warning('Invalid alpha convergence, using default value of %d.', alpha_conv); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% Verbosity
a = find(strcmp(varargin, 'verbose'));
verbose = true;
if(~isempty(a))
	v = varargin{a(1)+1};
	if(~isempty(v))
		verbose = logical(v);
	else
		warning('Invalid verbosity, using default value of ''true''.'); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% XScale
svals = {'linear', 'log'};
a = find(strcmp(varargin, 'xscale'));
xscale = 'linear';

if(~isempty(a))
	v = varargin{a(1)+1};
	if(~isempty(v))
		v = svals{strcmp(v, svals)};
	end
	
	if(~isempty(v))
		xscale = v;
	else
		warning('Invalid xscale, using default value of ''%s''', xscale); %#ok
	end
	
	varargin([a, a+1]) = [];		
end

% YScale
svals = {'linear', 'log'};
a = find(strcmp(varargin, 'yscale'));
yscale = 'linear';

if(~isempty(a))
	v = varargin{a(1)+1};
	if(~isempty(v))
		v = svals{strcmp(v, svals)};
	end
	
	if(~isempty(v))
		yscale = v;
	else
		warning('Invalid yscale, using default value of ''%s''', yscale); %#ok
	end
	
	varargin([a, a+1]) = [];		
end

% XLabel
a = find(strcmp(varargin, 'xlabel'));

if(~isempty(a))
	v = varargin{a(1)+1};
	
	if(ischar(v))
		xlab = v;
	end
	
	varargin([a, a+1]) = [];
end

% YLabel
a = find(strcmp(varargin, 'ylabel'));

if(~isempty(a))
	v = varargin{a(1)+1};
	
	if(ischar(v))
		ylab = v;
	end
	
	varargin([a, a+1]) = [];
end

% Data percentage
a = find(strcmp(varargin, 'dataperc'));
dataperc = 1e-4;
if(~isempty(a))
	dp = varargin{a(1)+1};
	if(~isempty(dp) && isnumeric(dp))
		dataperc = dp(1);
	else
		warning('Invalid data percentage, using default value of %d', dataperc); %#ok
	end
	
	varargin([a, a+1]) = [];	
end

% Optimization parameters
a = find(strcmp(varargin, 'optims'));
if(~isempty(a))
	os = varargin{a(1)+1};
	
	if(~isempty(os) && isstruct(os))
		o = os;
	else
		warning('Invalid optimset, usign default values.'); %#ok
	end
	
	varargin([a, a+1]) = [];
end

os = optimset('GradObj', 'on', 'Hessian', 'on',	'TolX', 1e-9, ...
	'TolFun',  1e-11, 'LargeScale', 'on',  'Display', 'off');

% Set default values where unset.
if(exist('o', 'var'))
	sn = fieldnames(os);
	
	nf = find(~isfield(o, sn));
	if(~isempty(nf))
		for i = 1:length(nf)
			o.(sn{i}) = os.(sn{i});
		end
	end
	
	for i = 1:length(sn)
		if(isempty(o.(sn{i})) && ~isempty(os.(sn{i})))
			o.(sn{i}) = os.(sn{i});
		end
	end
else
	o = os;
end

clear os;

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

% The "S" matrix is the diagonal matrix of singular values, organized in
% descending order, and are a representation of the magnitude of each
% eigenvalue in the Kernel vector. We can eliminate all combinations of S1
% and S2 which are below a certain fraction of the largest singular value
% (S1S2). These will be the singular values in our compressed space.
sd = diag(S1)*diag(S2)';
con = sd > S1(1, 1)*S2(1, 1)*dataperc;
n = sum(con(:));

S = sd(con);
[S, i] = sort(S, 'descend'); % Sort this like a proper SVD

% We want to vectorize our data matrix so that the calculations can go
% smoothly. We'll use the equality vec(K1*F*K2') = kron(K2, K1)*vec(F) to
% preserve the operations of the kernel in the vectorized space.
D = kron(D2, D1);
D = D(:, con);
D = D(:, i);   % Sort 

% U1 and U2 are the unitary transforms which transform between the SVD and
% data bases. Use an inverse transform on the data to bring the data into
% the SVD space, then select the subset of the data which are in our
% compressed space.
Mc = U1'*d.z*U2;
Mc = Mc(con);
Mc = Mc(i);		% MCI was a a telecommunications company before it changed 
					% its name to WorldCom. It is now owned by Verison.

% Create the kernel function in the compressed space by recombining S and
% D. S can be constructed this way because the singular value matrix is
% always a diagonal matrix with the singular values along the diagonal.
K = diag(S)*D';
id = eye(n, n);	% An identity matrix.

% Set up the output structure.
out.t1 = tau1;		
out.t2 = tau2;
out.f = {};
out.c = {};
out.K = K;
out.alpha = [];

ps = zeros(1, 2);
ps(1) = length(tau1);
ps(2) = length(tau2);

st = n*d.std;

% Now set up the guess vector.
a = find(strcmp(varargin, 'guess'), 1, 'first');
C0 = ones(n, 1); % Initial guess


if(~isempty(a))
	f = varargin{a+1};
	sf = size(f);

	
	if(isnumeric(f) && (all(sf == ps) || all(sf == fliplr(ps))))
		if(sf ~= ps)
			f = f';
		end
		
		C0 = -(K*reshape(f, prod(ps), 1) - Mc)/alphas(1);
		o.TypicalX = C0;
	else
		warning('The best guess provided was not valid, using default.'); %#ok
	end
end

ps = num2cell(ps);

% Subplot size
ss = ceil(sqrt(nAlphas));
conv = 0;

for i = 1:nAlphas
	out.alpha(i) = alpha;
	
	t = cputime;
	c = fminunc(@(c)fun(c, Mc, K, alpha, id), C0, o);
	
	t2 = cputime-t;
	minute = floor(t2/60);
	second = t2-minute*60;
		
	C0 = c;
	f = max(0, reshape(K'*c, ps{:}));
	
	out.c{i} = c;
	out.f{i} = f;
	
	if(verbose)
		fprintf('%d - Elapsed time %02.0g minutes %02.0g seconds.\n', i, minute, second);
	end
	
	subplot(ss, ss, i);
	contour(out.t1, out.t2, out.f{i}');

	set(gca, 'YScale', yscale);
	set(gca, 'XScale', xscale);

	title(sprintf('%d: \\alpha = %02.2g', i, alpha));
	
	if(exist('xlab', 'var'))
		xlabel(xlab);
	end
	
	if(exist('ylab', 'var'))
		ylabel(ylab);
	end
	
	drawnow;
	
	% Compute the new optimal alpha.
	a_opt = sqrt(st)/norm(c, 'fro');
	if(abs((alpha-a_opt)/alpha) < alpha_conv) % Convergence condition.
		conv = 1;
		break;
	end
	
	alpha = a_opt;
end

if(verbose)
	if(conv)
		fprintf('alpha converged at %3.3f after %d attempts.', alpha, length(out.alpha));
	else
		fprintf('alpha failed to converge after %d attempts.', nAlphas);
	end
end

if(conv)
	clf;
	
	contour(out.t1, out.t2, out.f{i}');

	set(gca, 'YScale', yscale);
	set(gca, 'XScale', xscale);

	title(sprintf('\\alpha = %02.2g', alpha));
	
	if(exist('xlab', 'var'))
		xlabel(xlab);
	end
	
	if(exist('ylab', 'var'))
		ylabel(ylab);
	end
end


function [f, G, H] = fun(c, Mc, K, alpha, id)
% Implements the unconstrained minimization from the referenced paper
% IEEE Transactions on Signal Processing, Vol. 50, No. 5 May 2002
%
% fr = max(0, K'c);

% This is the 'Kuhn-Tucker' condition.
% Here G(c) = K*min(0, Heavyside(diag(K'*c)))*K'
% The first two operators resolve to just K with the
% columns where K'*c <= 0 evalutes to true set to 0, so for speed
% we can use logical indexing and skip a matrix multiplication.
H = K;
H(:, K'*c <= 0) = 0;
H = H*K' + alpha*id;

G = H*c;	

f = 0.5*c'*G - c'*Mc;	% Eqn 31

if(nargout > 1)
	G = G - Mc;			% Eqn 32
end


