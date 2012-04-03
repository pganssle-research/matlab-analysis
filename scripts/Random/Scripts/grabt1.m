[num, txt, raw] = xlsread('C:\Users\Omega\Documents\T1Data.xlsx', 'DataGroup');

for i = 0:15
x(i+1, :) = num((i*8000+1):((i+1)*8000));
x(i+1, :) = x(i+1,:) - mean(x(i+1, :));
end

sr = 2000;
np = 8000;
x2 = zeros(8000, 1, 16);
x2(:, 1, :) = x(:, :)';

[c, times, spans] = get_mag(120, 1000/60, 50, sr, 25, x2, 0.85);

it = 0:(8/20):8;
c = -c;
[t1, t1_std] = find_t1_lsq(c', it')
calc_temp_from_T1(t1(1))