function o = fftx(varargin)

if(length(varargin) > 2)
    dim = varargin{3};
end

data = [varargin{1}];
s = size(data);

if(length(s) > 2)
    if(exist('dim', 'var'))
        n = s(dim);
    else
        n = dim(1);
    end
else
    n = length(data);
end
    
o = fft(data, varargin{2:end});
o = 2*o/n;
