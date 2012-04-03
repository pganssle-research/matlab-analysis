function hist = get_image_intensity_histogram(img, s1, s2) 
% Generates an intensity histogram for a given image (img) of size s1 x s2.
% This will simply be the integral within each square. This is essentially
% downsampling the image resolution.

s = size(img);

if(s(1) < s1)
    s1 = s(1);
    fprintf('Size 1 too large, setting to %d\n', s1);
end

if(mod(s(1), s1) ~= 0)
    s1 = s1-rem(s(1), s1);
    if(s1 <= 1)
        s1 = 1;
    end
    fprintf('Size 1 does not divide evenly into the size of the original image. Setting to %d\n', s1);
end

if(nargin < 3)
    s2 = 1;
end

if(s(2) < s2)
    s2 = s(2);
    fprintf('Size 2 too large, setting to %d\n', s2);
end

if(mod(s(2), s2 ~= 0))
    s2 = s2-rem(s(2), s2);
    if(s2 <= 1)
        s2 = 1;
    end
    fprintf('Size 2 does not divide evenly into the size of the original image. Setting to %d\n', s2);
end

% Break the image up into a cell containing the sub-arrays 
img = mat2cell(img, (s(1)/s1)*ones(s1, 1), (s(2)/s2)*ones(s2, 1));

% Generate the histogram.
hist = cellfun(@(x) mean(x(:)), img); 


