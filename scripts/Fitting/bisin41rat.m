function out = bisin41rat(x0, t)
% out = bisin41rat([w, A], t)
%
% out = A(sin(w*t)+sin(w*t/4)/2)

w = x0(1);
A = x0(2);

out = magnitude(A*(exp(-1i*w*t) + exp(-1i*w*t/4)))/2;

end