function index = outlier(data, alpha)
% Non-recursive outlier detection

if(~exist('alpha', 'var') || alpha < 0 || alpha > 1)
    alpha = 0.1; % 90% confidence level
end
y = data(:);

G = (mean(y)-y)/std(y);
n = length(y);

if(n < 3)
    index = zeros(length(y));
    return;
end

p = 1-(alpha/(2*n));
t = tinv(p, n-2);

index = G > ((n - 1)./(n.^(1/2)))*((t.^2)/(n-2+t.^2));
