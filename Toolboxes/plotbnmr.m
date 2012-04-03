function h = plotbnmr(A,varargin)
% PLOTBNMR	Plots NMR-spectra loaded from rbnmr-function
%
% SYNTAX	h = plotbnmr(A);
%
% IN		A:	Struct as read by rbnmr.
%			If A is not supplied, rbnmr is executed
%			in the current directory.
%
% OUT		h:	Handles to graphics created.
%			Each objekt in h is tagged with the
%			contents of an optional field A.Tag.

% Nils Nyberg, SLU, 2001-05-03
% nils.nyberg@kemi.slu.se
% rewritten, DFU, 2007-11-14
% nn@farma.ku.dk

%% Initialize some parameters
NoOfLevels = 6;		% No of plotting levels, if levels are calculated
SiNo = 4;			% Signal level over mean-signal, for calculating plotting levels


%% Check argument
if nargin == 0;
	disp('PLOTBNMR: Reads data from current directory');
	A = rbnmr;
end;

if nargin > 0;
	switch class(A);
		case 'cell';		% Loop and do recursive calls
			C = mybluered(length(A));
			H = cell(size(A));
			if ~ishold; hold on; end
			for i=1:length(A);
				if nargin == 1; 
					H{i} = plotbnmr(A{i},'-'); 
					set(H{i},'Color',C(i,:));
				else
					H{i} = plotbnmr(A{i},varargin{:});
				end
			end
			if nargout; h = cell2mat(H); end
			return;
		case 'struct';
			if ~all([isfield(A,'Data'), isfield(A,'XAxis')])
				error('PLOTBNMR: Argument has the wrong format');
			end
		otherwise
			error('PLOTBNMR: Argument must be a struct (or a cell of structs) from function rbnmr');
	end
end

%% 1D-plot
if ~isfield(A,'Proc2s')
	H = plot(A.XAxis,A.Data,varargin{:});
	set(gca,'XDir','reverse');
else

%% 2D-plot
if ~isfield(A,'Levels')
    A.Levels = A.Procs.S_DEV*SiNo*[1 cumprod(ones(1,NoOfLevels-1)*1.1)];
end

if (length(varargin) == 1);
    [C,H] = contour3(A.XAxis,A.YAxis-mean(diff(A.YAxis))/2,A.Data,A.Levels,varargin{:});
else
    [C,H] = contour3(A.XAxis,A.YAxis-mean(diff(A.YAxis))/2,A.Data,A.Levels);
end;

view(2);grid off;box on;
set(gca,'XDir','reverse','YDir','reverse','YAxisLocation','right');
end


%% Finish
if isfield(A,'Tag');
	set(H,'Tag',A.Tag);
end;
try
	set(H,'DisplayName',A.Info.PlotLabel);
catch ME
	try
		set(H,'DisplayName',A.Title);
	catch ME2
		% Do nothing
	end
end

if nargout == 1; 
	h = H;
	if numel(h) > 1
		H = hggroup;
		set(h,'Parent',H);
	end
end;

function c=mybluered(n)
% My blue to red colormap
n = ceil(n);
if nargin < 1; n = 64; end
if mod(n,2); n = n + 1; end

hsv = [...
	[repmat(0.6,n/2,1) linspace(1.0,0.6,n/2)' repmat(1.0,n/2,1)];...
	[repmat(1.0,n/2,1) linspace(0.6,1.0,n/2)' repmat(1.0,n/2,1)]];

c = hsv2rgb(hsv);





