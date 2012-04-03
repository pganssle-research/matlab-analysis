function Dw = find_diff(curve, V, G, t, tau, temp)

if(~exist('temp', 'var'))
    temp = 60;
end

Dw_est = calc_Dw_from_temp(temp);

sc = size(squeeze(curve));
sv = size(squeeze(V));

% I want these to both be vectors of size n x 1
if(length(sc) > 2) 
   error('Invalid size of curve. Must be n x 1');
end

if(length(sv) > 2) 
    error('Invalid size of voltage vector. Must be n x 1');
end

if(sc(2) > 1)
    curve = curve';
    sc = size(curve);
end

if(sc(2) > 1)
   error('Must be 1D data.'); 
end

if(sv(2) > 1)
    V = V';
    sv = size(V);
end

if(sv(2) > 1)
    error('Must be 1D voltage vector.');
end

if(sv(1) ~= sc(1))
    error('Vectors must be the same size.');
end

o = optimset('Display', 'off');  
S = -(2*pi*4257 * G)^2 * tau^2 * t/3;

Dw = lsqcurvefit(@(x, t)x(2)*exp(S*x(1)*t.^2), [Dw_est, curve(1)], V, curve, [], [], o);