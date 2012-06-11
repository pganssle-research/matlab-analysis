function out = read_vnmr_data(filename, lb)
% Reads in varian data. Code is mostly stolen from matnmr.
%
%
% Returns a struct with the format:
%
% out.fid => The FID
% out.t   => Time vector
% out.fft => fourier transform, size is the next highest power of 2.
% out.w   => Frequency vector (Hz)
% out.p   => Frequency vector (ppm)
% out.np  => Number of points
% out.sw  => Spectral width
% out.at  => Acquisition time
% out.sfrq=> Spectrometer frequency
%
% In case of a blank filename, use a folder select dialog
%
% out = read_vnmr_data(filename, lb);


% Get the directory containing the filename, if one is not supplied
if(nargin < 1)
    if(exist('uigetdir2.m', 'file'))
        filename = uigetdir2;
    else
        filename = uigetdir;
    end
end

% Make sure it has a fileseparator as a terminating character.
if(filename(end) ~= filesep)
    filename(end+1) = filesep;
end

% Get our filenames to be read in
ppfile = [filename 'procpar'];
fidfile = [filename 'fid'];

% Check that this is a valid file.
if(~exist(ppfile, 'file'))
    error('Procpar file not found');
elseif (~exist(fidfile, 'file'))
    error('FID File not found.');
end

%Read in the parameters from the procpar file.
pp_read = ReadParameterFile(ppfile);  % Reads in the full procpar file

% Grab some parameters.
pn = 'pp_read';     % Convenience
eval(gen_extraction_func('sfrq', pn, 'out.sfrq')); % Spectrometer freq.
eval(gen_extraction_func('sw', pn, 'out.sw')); % Spectral width
eval(gen_extraction_func('np', pn, 'out.np')); % Number of points
eval(gen_extraction_func('at', pn, 'out.at')); % Acquisition time.

clear pp_read;

% Read the file itself.
fp = fopen(fidfile, 'r', 'b'); % Big-endian

% First read the header out into a struct
h = struct('nblocks', 0, 'ntraces', 0, 'np', 0, 'ebytes', 0, 'tbytes', 0, 'bbytes', 0, 'vers_id', 0, 'stat', 0, 'nbheaders', 0);
[a, ~] = fread(fp, 6, 'int32');
h.nblocks = a(1);
h.ntraces = a(2);
h.np = a(3);
h.ebytes = a(4);
h.tbytes = a(5);
h.bbytes = a(6);

[a, ~] = fread(fp, 2, 'int16');
h.vers_id = a(1);
h.stat = a(2);

[a, ~] = fread(fp, 1, 'int32');
h.nbheaders = a(1);

% Status, bitwise breakdown: (VNMR Progrogramming manual)
% Changed from a 0-based index to reflect Matlab's 1-based index.
% Bits 1-7: File header and block header status bits (7 is unused).
% 1 S_DATA          0x1 0 = no data, 1 = data 
% 2 S_SPEC          0x2 0 = FID, 1 = spectrum
% 3 S_32            0x4 *
% 4 S_FLOAT         0x8 0 = integer, 1 = floating point
% 5 S_COMPLEX       0x10 0 = real, 1 = complex
% 6 S_HYPERCOMPLEX  0x20 1 = hypercomplex
%
% Bits 8-16: File header status bits (11 and 16 are unused)
% 8 S_ACQPAR        0x80 0 = not Acqpar, 1 = Acqpar 
% 9 S_SECND         0x100 0 = first FT, 1 = second FT
% 10 S_TRANSF        0x200 0 = regular, 1 = transposed
% 12 S_NP           0x800 1 = np dimension is active
% 13 S_NF           0x1000 1 = nf dimension is active
% 14 S_NI           0x2000 1 = ni dimension is active
% 15 S_NI2          0x4000 1 = ni2 dimension is activ

% Use this to determine the binary format

if(~bitget(h.stat, 3) && ~bitget(h.stat, 4))
    prec = 'int16';
elseif(bitget(h.stat, 3) && ~bitget(h.stat, 4))
    prec = 'int32';
else
    prec = 'float';
end

% Now read the data.
np = h.np;
if(np ~= out.np)
    out.np = np;
    disp('Warning: Number of points in procpar does not match binary header.');
end

% Pre-allocate a matrix.
data = zeros(h.nblocks*h.ntraces*np, 1);

for n = 1:h.nblocks
    % Read in the headers (and throw them away)
    for m = 1:h.nbheaders
        fread(fp, h.nbheaders*14, 'int16'); % If nbheaders >= 1, this might double-count...
    end

    % Grab the actual data.
    [tmp, ~] = fread(fp, h.ntraces*np, prec);

    % Zero pad if the length is wrong for whatever reason.
    if(length(tmp) < h.ntraces*np)
        tmp((length(tmp)+1):h.ntraces*np) = 0;
        disp('Didn''t read enough values from file, zero-packing.');
    end

    data(((n-1)*h.ntraces*np+1):(n*h.ntraces*np)) = tmp;
end

fclose(fp);

% Put the data into a useful form. This involves reshaping this into a 2D
% array, where the number of blocks (?) and traces (presumably transients?)
% are along one dimension and the FIDs are along the other. These are
% complex spectra, stored in a striped manner, so each 2-bit pair in the
% data array is one complex number.

 out.fid = zeros(np/2, h.nblocks*h.ntraces);
 data = reshape(data, np, h.nblocks*h.ntraces);
 out.fid(1:(np/2), :) = data(1:2:np, :) + 1i*data(2:2:np, :);
 
 % Update the number of points.
 out.np = np/2;
 np = np/2;
 
 % Generate the other vectors
 out.t = ((0:(np-1))/out.sw)'; % Time vector (seconds)
 
 % Add apodization to structure
 if(exist('lb', 'var') && lb > 0)
	out.lb = lb;
else
	out.lb = 0;
end

out = mouse_fft(out, 1);
 
function funs = gen_extraction_func(pname, param_ar, func_out)
% This function is a bit bizarre. Since Matlab is pass-by-value, I am
% hesitant to pass the full parameter file (which will be a large cell
% array). However, it is convenient to have a centralized way of doing
% this, so if you want to get the value for 'sw' in parameter cell array
% 'pp_read', and put it into variable, 'out.sw', you would use:
%
% eval(gen_extraction_func('sw', 'pp_read', 'out.sw'));
%
% Depreciated because it was a dumb idea in the first place. 

l = cell(11, 1);
nl = sprintf('\n');

l{01} = ['cout = '''';' nl];
l{02} = ['fstr = ''' pname ' '';' nl];
l{03} = ['loc = find(strncmp(fstr, ' param_ar ', length(fstr)));' nl];
l{04} = ['if(~isempty(loc))' nl];    % If it's empty, we'll return empty
l{05} = [    'nstr = deblank(' param_ar '(loc+1));' nl];
l{06} = [    'vals = sscanf(char(nstr), ''%d %lf'');' nl];
l{07} = [    'if(length(vals) >= 2)' nl];
l{08} = [        'cout = vals(2);' nl];
l{09} = [    'end' nl];
l{10} = ['end' nl];

if(nargin > 2)
    l{11} = [func_out ' = cout;'];   % This is how we'll return the value.
else
    l{11} = ['cout;'];
end

% Generate the function (probably can be vectorized better, but char(l)
% does not work well.).
funs = '';
for i = 1:11
    funs = sprintf('%s%s', funs, l{i});
end

function Parameters = ReadParameterFile(FileName)
% Generates a cell array of all the lines in the procpar file.

fp = fopen(FileName, 'r');
if (fp == -1)
  disp(['Couldn''t open parameter file "' FileName '". Aborting ...']);
  Parameters = [];
  return
end

%
%read first line and edit
%
tmp = fgetl(fp);
tmp([strfind(tmp, '#') strfind(tmp, '$') strfind(tmp, ';')]) = '';
Parameters = {tmp};

while 1
  tmp = fgetl(fp);
  if ischar(tmp)
    tmp = deblank(tmp);
    if (length(tmp)<80)		%avoid endless lines and thus and memory time loss
      tmp([strfind(tmp, '#') strfind(tmp, '$')]) = '';
      Parameters = [Parameters, {tmp}];
    end
  else
    break
  end
end

fclose(fp);











