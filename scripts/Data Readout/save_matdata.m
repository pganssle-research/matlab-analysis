function save_matdata(struct, path)
% Function for saving data as a .mat for export to pylab.
%
% Either you can either pass a full path or a name to 'path'
%
% save_matdata(struct, path)

s = struct;

% Parse the input arguments
if(exist('path', 'var'))
	% Assume it's a name.
	[p,n,e] = fileparts(path);
	
	if(isempty(e))
		if(~isempty(p))
			p = [path filesep n];
			n = '';
		else
			n = [n '.mat'];
		end
	else
		n = [n e];
	end
	
	if(exist(p, 'file'))
		path = p;
	end
else
	path = '';
	n = '';
end

if(isempty(path))
	path = get_path('save_path.mat', '.mat', 0, n);
end

t = s.t;
md = s.mdata;
ad = s.adata;

f = '';
if(isfield(s, 'f'))
	f = s.f;
end

sp = [];
if(isfield(s, 'fft'))
	sp = s.fft;
end

% Calibration to magnetic field.
mag_cal = 1.0;
if(isfield(s, 'disp') && isfield(s.disp, 'mag_cal'))
	mag_cal = s.disp.mag_cal;
end

% Check for fits
t1 = [];
t1e = [];
t1a = [];
t1ae = [];

t2 = [];
t2e = [];
t2a = [];
t2ae = [];

it = [];

V = [];
G = [];
M = [];

D = [];
De = [];
n = [];
tau = [];

std = [];
c = [];
co = [];
cf = [];

if(isfield(s, 'fit'))
	% Posible field names:
	{'t1', 't1e', 't1a', 't1ae', 't2', 't2e', 't2a', 't2ae', 'V', ...
		'D', 'De', 'Da', 'Dae', 'n', 'tau', 't', 'c', 'co', 'cf', 'c'};


	
	
	
	
	
	
	