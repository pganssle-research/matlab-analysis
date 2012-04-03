function f = polyphase(p, freq)
% Returns an exponential e^(i*polyval(p, freq))
%
% Usage:
% o = polyphase(p, freq);

f = exp(1i*polyval(p, freq));
   