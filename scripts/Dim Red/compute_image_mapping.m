function [mappedA, mapping] = compute_image_mapping(img, type, varargin)
% Takes images and vectorizes them appropriately to be passed to
% compute_mapping.

img = img(:, :); % Easy vectorization
[mappedA, mapping] = compute_mapping(img, type, varargin{:});