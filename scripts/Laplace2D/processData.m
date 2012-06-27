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
% 'alphaLims' -> Limits (in logspace) to try. Default: [-5, 5]
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
%					'Algorithm' = 'interior-point'
%					'GradObj' = 'on'
%					'Hessian' = 'off'
%					'TolX' = 1e-9
%					'TolFun' = 1e-12
%					'MaxIter' = 5000
%					'MaxFunEvals' = 5000
%					'Display', 'off'
%
%					The minimization function can return Hessians, though this
%					is not supported by the interior-point or sqp algorithms.
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

% Alpha limits
a = find(strcmp(varargin, 'alphaLims'), 1, 'first');
alphaLims = [-10, 10];
if(~isempty(a))
	al = varargin{a(1)+1};
	if(isnumeric(al) && length(al) == 2)
		alphaLims(:) = al(1:2);
	else
		warning('Invalid alphaLims, using default value of [%d, %d].', alphaLims(1), alphaLims(2)); %#ok
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

os = optimset('Algorithm', 'interior-point', 'GradObj', 'on', 'Hessian', ...
	'off', 'MaxIter', 5000, 'TolX', 1e-8, 'LargeScale', 'on', ...
	'TolFun',  1e-12, 'MaxFunEvals', 5000, 'Display', 'off');

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
[S, i] = sort(S, 'descend');

% We want to vectorize our data matrix so that the calculations can go
% smoothly. We'll use the equality vec(K1*F*K2') = kron(K2, K1)*vec(F) to
% preserve the operations of the kernel in the vectorized space.
D = kron(D2, D1);
D = D(:, con);
D = D(:, i);

% U1 and U2 are the unitary transforms which transform between the SVD and
% data bases. Use an inverse transform on the data to bring the data into
% the SVD space, then select the subset of the data which are in our
% compressed space.
Mc = U1'*d.z*U2;
Mc = Mc(con);
Mc = Mc(i);

% Create the kernel function in the compressed space by recombining S and
% D. S can be constructed this way because the singular value matrix is
% always a diagonal matrix with the singular values along the diagonal.
K = diag(S)*D';
id = eye(n, n);	% An identity matrix.
Ms = Mc'*Mc;		% Precalculate the norm of the compressed data for speed

% Set up the output structure.
out.t1 = tau1;		
out.t2 = tau2;
out.f = {};
out.F = {};
out.K = K;

% Now set up the guess vector.
a = find(strcmp(varargin, 'guess'), 1, 'first');
F0 = ones(n, 1); % Initial guess

ps = zeros(1, 2);
ps(1) = length(tau1);
ps(2) = length(tau2);

if(~isempty(a))
	f = varargin{a+1};
	sf = size(f);

	
	if(isnumeric(f) && (all(sf == ps) || all(sf == fliplr(ps))))
		if(sf ~= ps)
			f = f';
		end
		
		F0 = K*reshape(f, prod(ps), 1);
		o.TypicalX = F0;
	else
		warning('The best guess provided was not valid, using default.'); %#ok
	end
end
ps = num2cell(ps);

alphas = logspace(alphaLims(1), alphaLims(2), nAlphas);
out.alphas = alphas;

% Subplot size
ss = ceil(sqrt(nAlphas));
H1 = K*K';

for i = 1:nAlphas
	alpha = alphas(i);
	H = alpha*H1;
	
	t = cputime;
	
	F = fmincon(@(F)fun(F, Mc, H, id, Ms), F0, [], [], [], [], zeros(size(F0)), [], [], o);
	
	t2 = cputime-t;
	minute = floor(t2/60);
	second = t2-minute*60;
	fprintf('Elapsed time %02.0g minutes %02.0g seconds.\n', minute, second);
	
	out.F{i} = F;
	out.f{i} = reshape(K'*F, ps{:});
	
	subplot(ss, ss, i);
	image(out.t1, out.t2, out.f{i}', 'CDataMapping', 'scaled');
	set(gca, 'YDir', 'normal');
	
	title(sprintf('\\alpha = %d', alpha));
	
	xlabel('T1');
	ylabel('T2');
	
	drawnow;
end

function [f, g, H] = fun(F, Mc, H, id, Ms)
% This calculates the minimization function:
%
% f = ||Mc - F||^2 + alpha*||K'*F||^2
%   = sqrt(sum(Mc.^2 + F.^2 - 2*F.*Mc)).^2 + alpha*sqrt(sum(K'*F))^2
%   = abs(Mc'*Mc + F'*(F-Mc)) + alpha*F'*K*K'*F;
%   = abs(Ms + F'*(F-Mc)) + F'*H*F;
%
% g = Grad = df/dF
%   = 2*(F - Mc + alpha*K*K'*F);
%   = 2*(F - Mc + H*F)
%
% H = Hess = dg/dF = 2*(alpha*K*K' + id) = 2*(H + id);  
%
% Where F is the spectrum in the compressed data set.
%
% This is the most speed-critical part of the algorithm, as it is
% executed as many as MaxIter times.

G1 = F - 2*Mc;
G2 = H*F;

f = abs(F'*G1 + Ms) + F'*G2;

if(nargout > 1)
	g = G1 + 2*G2 + F;
end

if(nargout > 2)
	H = 2*(H+id);
end

