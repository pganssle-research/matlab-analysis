function [out, stdev] = make_avg(d, t_vec)
% Input curves appropriately.

data = d;

s = size(data);
new_data = zeros(s);
if nargin == 1
    t_vec = transpose(1:s(1));
end

o = optimset('Display', 'off');
t1_base = lsqcurvefit(@exponential_fit, [5.0, 1.0], t_vec, data(:, 1), [], [], o);
t1_base = t1_base(2);

for i=2:s(2)
    t1 = lsqcurvefit(@exponential_fit, [5.0, 1.0], t_vec, data(:, i), [], [], o);
    t1 = t1(2);
    
    data(:, i) = data(:, i)*t1_base/t1;
end


out = zeros(s(1), 1);
stdev = zeros(s(1), 1);

for i=1:s(1)
    out(i) = mean(data(i, :));
    stdev(i) = std(data(i, :));
end

