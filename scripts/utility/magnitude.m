function out = magnitude(in)
% Insert a complex vector, it gives you out the magnitude.
%
% out = (real(in).^2 + imag(in).^2).^(1/2);
%
% Usage:
% out = magnitude(in);

out = (real(in).^2 + imag(in).^2).^(1/2);