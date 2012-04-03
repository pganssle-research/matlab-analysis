function [f, s, sa] = get_spectrum(data, po, cut, apo)
% Gets the poly-subtracted, apodized and cut-and-zero-packed version
% spectrum from a TDMExperiment structure.
%
%
% Usage:
% get_spectrum(data, po, cut, apo);

if(isstruct('data'))
    error('Data must be structure.');
end

sr = data.prog.sr;
data = data.mdata;
si = num2cell(size(data));

if(~exist('cut', 'var'))
    cut = -1;
end

if(~exist('po', 'var'))
   po = -1; 
end

if(~exist('apo', 'var'))
    apo = -1;
end

if(po >= 0)
    data = sub_poly(data, po, 0, cut);
end

if(cut > 0)
    data = cutandzeropack(data, cut);
end

if(apo > 0)
   data = apodize1d(data, apo, sr); 
end

nnp = 2^(ceil(log2(size(data, 1))));
f = linspace(0, sr/2, nnp/2);
s = fft(data, nnp);
s = s(1:(nnp/2), :);
s = reshape(s, nnp/2, si{2:end});

sa = mean(s, 2);

if(length(si) > 2)
    sa = reshape(sa, nnp/2, si{3:end});
end