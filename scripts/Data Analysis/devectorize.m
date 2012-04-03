function devec = devectorize(vec, os)
% Reverse operation of "vectorize" - reshapes a vectorized image to the
% original size. Feed this the outputs of vectorize.
%
% devec = devectorize(vec, os);

devec = reshape(vec, os);