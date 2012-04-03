samps;

s = s/100;
n = length(s);
img = zeros(n, 128, 384);

g1 = gauss2d(5, [0, 64], 128, 384, -64);
g2 = gauss2d(5, [0, 64+128], 128, 384, -64);
g3 = gauss2d(5, [0, 64+256], 128, 384, -64);

x = s(:, 1);
y = s(:, 2);
z = s(:, 3);

if(~exist('glu_s', 'var'))
    for i = 1:n
        img(i, :, :) = g1*x(i) + g2*y(i) + g3*z(i);
    end
    img = reshape(img, n, 128*384);
else
   img = zeros(n, length(glu_s));
   for i = 1:n
       img(i, :) = glu_s*x(i) + lys_s*y(i) + pro_s*z(i);
   end
   
   ref = arrayfun(@(j)logical(j > 1 && j < 4.5), p);
   img = img(:, ref);
end

figure(2)
scatter3(x, y, z);
axis equal

figure(1)


% Center data and set unit variance.
img = bsxfun(@minus, img, mean(img));    % Zero mean
img = bsxfun(@rdivide, img, (std(img)~=0).*std(img) + (std(img)==0));   % Unit variance

dims = 2;

[pc_a, pmap] = compute_mapping(img, 'PCA', dims);

if(~exist('cols', 'var'))
    cols = 1:length(x);
end

if(dims > 2)
    scatter3(pc_a(:, 1), pc_a(:, 2), pc_a(:, 3), cols);
else
    scatter(pc_a(:, 1), pc_a(:, 2), 50, cols, 'filled');
end
axis equal