function o = mc_options(varargin)
% Generates an options set for use with mc_struct_viewer. Any arguments
% not found will be set to default.
%
% Options:
% .fft:		0 - > Start with FID (default)
%				1 - > Start with real channel
%				2 - > Start with imaginary channel.
%				3 - > Start with the magnitude
%
% c:			0 - > View the FID/FFT (default)
%				1 - > View the curves array (lobe1-lobe2, lobe2-lobe3, etc)
%				2 - > View the points array (lobe1, lobe2, lobe3, etc.)
%
% step:		 0 - > Show every transient (equivalent to (:)).
%				-1 - > Show the average (default)
%				-2 - > Show the first one
%
% newfig:	Bool - > Whether or not to open in a new figure. (Default: 0)
% name:		Str - > Name of the figure (Default: 'MC Data Viewer') 
%
%
% Usage:
% o = mc_options('param1', 'val1', ...);

o = struct('fft', 0, 'c', 0, 'step', [-1], 'newfig', 0, 'name', 'MC Data Viewer');

% Types:
TYPE_NUM = 0;
TYPE_BOOL = 1;
TYPE_STR = 2;

options = {'fft', 'c', 'step', 'newfig', 'name'};
type = {TYPE_NUM, TYPE_NUM, TYPE_NUM, TYPE_BOOL, TYPE_STR};

maxes = {3, 3, [], [], []};
mins = {0, 0, -2, [], []};

n = floor(length(varargin)/2);

for i = 1:length(options)
	if(n < 1)
		break;
	end
	
	ind = find(strcmp(options(i), {varargin{1:2:end}}), 1, 'first');
	ind = ind*2 - 1; %Odds are always parameters.
	if(~isempty(ind))
		val = varargin{ind+1};
		
		if(type{i} == TYPE_NUM && isnumeric(val))
			if((isempty(maxes{i}) || val <= maxes{i}) && (isempty(mins{i}) || val >= mins{i}))
				o.(varargin{ind}) = val;
				n = n-1;
			end
		elseif(type{i} == TYPE_BOOL)
			o.(varargin{ind}) = logical(val);
			n = n-1;
		elseif(type{i} == TYPE_STR && ischar(val))
			o.(varargin{ind}) = val;
			n = n-1;
		end
	end
end
