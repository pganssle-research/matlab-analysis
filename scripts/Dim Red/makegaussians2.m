n = size(r, 2);
img = zeros(n, 128, 384);

% gaussian = gauss2d(sigma, center_point, size_1, size_2, [bottom, left, spacing_1, spacing_2])
% In this case, it generates matrices of size 128 x 384 which go from 0->384 and -64->64
g1 = gauss2d(5, [0, 64], 128, 384, -64);
g2 = gauss2d(5, [0, 64+128], 128, 384, -64);
g3 = gauss2d(5, [0, 64+256], 128, 384, -64);

% g1 = g1/sum(sum(g1));
% g2 = g2/sum(sum(g2));
% g3 = g3/sum(sum(g3));

% g1 = g1/max(max(g1));
% g2 = g2/max(max(g2));
% g3 = g3/max(max(g3));

% for i = 1:7
% img(i, :, :) = bitget(i, 1)*g(1, :, :) + bitget(i, 2)*g(2, :, :) + bitget(i, 3)*g(3, :, :);
% end
   
% k = 0;
% for t = linspace(0, 2*pi, nt)
%     for p = linspace(-pi/2, pi/2, np)
%         k = k+1;
%         img(k, :, :) = g1*sin(t)*cos(p) + g2*sin(t)*cos(p) + g3*cos(t);
%     end
% end

% Spiral method
% phi = (1+sqrt(5))/2;
% angle = 0;
% x = zeros(1, n);
% y = zeros(1, n);
% z = linspace(-1, 1, n);
% for i = 1:n
%     r = sqrt(1-z(i)^2);
%     x(i) = cos(angle)*r;
%     y(i) = sin(angle)*r;
%     
%     img(i, :, :) = g1*x(i) + g2*z(i) + g3*y(i);
%     angle = angle+phi;
% end

% Uniform random distribution
% xyz = (rand(n, 3)*2)-1;
% xyz = (xyz./(repmat(arrayfun(@(x)norm(xyz(x, :)), 1:size(xyz, 1)), 3, 1))');
% x(:) = xyz(:, 1);
% y(:) = xyz(:, 2);
% z(:) = xyz(:, 3);

% Uniform random distribution 2
% [x,y,z] = sph2cart(pi*rand(n, 1), 2*pi*rand(n, 1), ones(n, 1));

n = size(r, 2);
x = r(1, :);
y = r(2, :);
z = r(3, :);

for i = 1:n
    img(i, :, :) = g1*x(i) + g2*y(i) + g3*z(i);
end

figure(2)
colormap(hsv);
scatter3(x, y, z, 30, z);
axis equal

figure(1)

img = reshape(img, n, 128*384);

% Center data and set unit variance.
img = bsxfun(@minus, img, mean(img));    % Zero mean
img = bsxfun(@rdivide, img, (std(img)~=0).*std(img) + (std(img)==0));   % Unit variance

[pca, pmap] = compute_mapping(img, 'PCA', 2);
colormap(hsv)

dim = size(pca, 2);

if(dim == 1)
    scatter(pca(:), ones(size(pca)), 30, z);
elseif (dim == 2)
   scatter(pca(:, 1), pca(:, 2), 30, z); 
else
   scatter3(pca(:, 1), pca(:, 2), pca(:, 3), 30, z);
end

axis equal