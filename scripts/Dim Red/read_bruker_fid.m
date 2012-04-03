function data = read_bruker_fid(path)

if(~exist('path', 'var') || ~exist(path,'file'))
    path = uigetdir2('C:\Users\Omega\Dropbox\');
end

[names, paths] = getAllFiles(path);

proc_ind = find(strcmp(names, 'procs'), 1, 'first');
if(isempty(proc_ind))
    error('Procpar file missing.')
end

ProcPar = ReadParameterFile2(paths{proc_ind});

% Get the number of points in the main dimension
si = 0;
[buff, count] = get_proc_param(ProcPar, 'SI', 'd');

if(count == 0)
    error('No points indicated in the procpar. Can''t proceed.');
else
    si = buff;
end

% Get the spectral width
sw = 0.0;
[buff, count] = get_proc_param(ProcPar, 'SW_p', 'lf');

if(count == 0)
    error('No spectral width could be found. Can''t proceed.');
else
    sw = buff;
end

% Get the byte order
big_endian = logical(0);
[buff, count] = get_proc_param(ProcPar, 'BYTORDP', 'd');

if(count == 0)
    warning('Byte order missing, assuming little endian. Wish us luck.')
else
   big_endian = logical(buff); 
end

% Get the spectrometer frequency
sf = 400.0;
[buff, count] = get_proc_param(ProcPar, 'SF', 'f');

if(count == 0)
    warning('Spectrometer frequency missing. Assuming 400.0 MHz')
else
   sf = buff;
end

% Get the dimensionality (Currently unused)
dim = 1;
[buff, count] = get_proc_param(ProcPar, 'PPARMOD', 'd');

if(count ~= 0)
    dim = buff+1;
end

% Read the FID file
fid_ind = find(strcmp(names, 'fid'), 1, 'first');
if(isempty(fid_ind))
    error('FID file missing.')
end

fp = fopen(paths{fid_ind});

if(~big_endian)
    mformat = 'b';
else
    mformat = 'l';
end

fid_data = fread(fp, si, 'int32', 0, mformat);
fclose(fp);

fid_data = fid_data';
t = (0:(si-1))/sw;

n_ind = fliplr(find(path == '\', 2, 'last'));

if(n_ind(1) < length(path))
   title = path((n_ind(1)+1):end);
else
    title = path((n_ind(2)+1):end);
end

data.name = title;
data.fid = fid_data;
data.t = t;

% Do some basic data processing
m_t = max(data.t);
np_fft = 2^(ceil(log2(si))+1);

% Shift everything so that the first point is the one with the largest
% absolute value.
[~, i] = max(abs(fid_data));
fid_data = fid_data(i:end);
t = t(i:end);

fid_data = fid_data .* exp(-3 * t/m_t); % Apodize the data

fft_data = fft(fid_data, np_fft);
f = linspace(0, sw, np_fft);
w = f + sf;                     % In absolute frequency units
p = w/sf;

data.f = f;
data.w = w;

data.p = p;
data.fft = fft_data;

data.sf = sf;
data.sw = sw;

% Processing from the procpar

% Frequency offset
[buff, count] = get_proc_param(ProcPar, 'OFFSET', 'f');
if(count ~= 0)
    data.p = data.p - buff;
end

% Phase corrections
[buff, count] = get_proc_param(ProcPar, 'PHC0', 'f');
phc0 = 0.0;
if(count ~= 0)
    phc0 = buff;
end

[buff, count] = get_proc_param(ProcPar, 'PHC1', 'f'); 
phc1 = 0.0;
if(count ~= 0)
    phc1 = buff;
end

function [buff, count] = get_proc_param(ProcPar, param, type)
% Type is d or f.
buff = [];
count = 0;
buff_ind = find(strncmp(ProcPar, [param, '='], length(param)), 1, 'first');

if(~isempty(buff_ind))
    [buff, count] = sscanf(ProcPar{buff_ind}, ['%*s%' type]);
end



