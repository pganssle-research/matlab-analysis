function [out, spans, subset, tc] = get_subset(data, len, start, asym, frac, num_windows, sr)
% Gets a subset of the data, returns a variable "spans" indicating where
% the subsets were taken from (scaled to the min/max of the data).
%
% All outputs other than data and len are optional. sr is only optional if
% data is a standard struct. Pass 'adata' to 'sr' if you want to use the
% average data from the struct.
%
% start, asym, frac and num_windows will use default values if they are set
% to any negative number.
%
% Default values:
% start = 0 ms
% asym = 0.5
% frac = 0.5
% num_windows = (all available)
%
% Usage:
% [out, spans, subset, tc] = get_subset(data, len, start, asym, frac, num_windows, sr);

if(~exist('data', 'var') || (isstruct(data) && ~isfield(data, 'mdata') && ~isfield(data, 'adata')))
    error('Must supply data!');
end

if(~exist('len', 'var'))
    error('Must supply a length of the subset.');
end

adata = 0; % Boolean, for later

if(~exist('sr', 'var') || ischar(sr))
    if(exist('sr', 'var') && strcmp(sr, 'adata'))
        adata = 1;
    end
    
    if(~isstruct(data) || ~isfield(data, 'prog') || ~isstruct(data.prog) || ~isfield(data.prog, 'sr'))
        error('Must provide sampling rate if data is not a well-formed structure.');
    else
       sr = data.prog.sr; 
    end
end

sr = sr/1000;  % Convert to kHz.

% We need to read the data we'll be working with for the next part
if(isstruct(data))    
   if(adata && isfield(data, 'adata'))
       dat = data.adata;
   else
       dat = data.mdata;
   end
else
    dat = data;
end

dat = squeeze(dat);

if(isempty(dat))
    error('Must supply data!');
end

% Data can be any size, but it must be of the form [data {everything else}].
ds = size(dat);
tl = ds(1)/sr;     % Total length, in ms.

if(len > tl)
   error('Subset length cannot be longer than total length.'); 
end

if(~exist('start', 'var') || start < 0 || start >= tl)
    start = 0;
end

tl2 = tl-start;     % Total length, minus start offset

start = start*sr; % How many samples is this?

if(~exist('num_windows', 'var') || num_windows > tl2/len || num_windows <= 0)
    % If num_windows is not specified or is too many, get them all.
    num_windows = floor(tl2/len);   
end

dlen = len*sr;
if(dlen ~= round(dlen))
 %  warning('Length is not an even multiple of the sampling rate - this could cause minor issues.'); 
end

dlen = len*sr; 
if(dlen <= 0)
    error('Length too short.');
end

if(~exist('frac', 'var') || frac <= 0 || frac > 1)
    frac = 0.5; % Make it about half the length.
end

window = ceil(frac*dlen); % How many points in each window.

if(window == 0)
    error('Fraction of length too short.');
end

if(~exist('asym', 'var') || asym <= 0 || asym > 1)
   asym = 0.65; 
end

% Now we have the inputs we need, we can get the outputs.
if(length(ds) > 1)
    c = num2cell(ds(2:end));
else
    c = {1};
end

tlen = length(dat(1, :));
spans = zeros(ds(1), 1);

% Where to start within a window
off = (dlen-window)*asym;
if((off+window)>dlen)
    off = dlen-window;
end

% Where to sample from
indices_pos = cell2mat(arrayfun(@(x)(x*dlen+off+1):(x*dlen+off+window), 0:2:(num_windows-1), 'UniformOutput', false))' + start;
indices_neg = cell2mat(arrayfun(@(x)(x*dlen+off+1):(x*dlen+off+window), 1:2:(num_windows-1), 'UniformOutput', false))' + start;
indices = cell2mat(arrayfun(@(x)(x*dlen+off+1):(x*dlen+off+window), 0:1:(num_windows-1), 'UniformOutput', false))' + start;

% Round at the end to reduce rounding errors.
indices_pos = round(indices_pos);
indices_neg = round(indices_neg);
indices = round(indices);

% Set the span outputs - scaled to the outputs
% Now the actual output vector
out = cell2mat(arrayfun(@(x)dat(indices, x), 1:tlen, 'UniformOutput', false));
out = reshape(out, window, num_windows, c{:});

Min = min(out(:));
Max = max(out(:));
s = (Max-Min)*0.05;
smean = mean(out(:));

spans(:) = smean;
spans(indices_neg) = Min-s;
spans(indices_pos) = Max-s;

subset = indices;
tc = (indices(1:window:(end-1))+indices(window:window:end))/(2*sr);






