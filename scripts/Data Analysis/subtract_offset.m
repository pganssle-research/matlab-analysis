function out = subtract_offset(in)
% Function that subtracts a DC offset from your data.
% 
% Usage: out = subtract_offset(in);
%
% This assumes that this is for a 2D matrix, pointsxtransients
% To be updated later for multiple dimensions.

s = size(in);
nd = 1;
if(length(s) > 2)
    nd = length(s)-1;
end
out = zeros(s);

for i = 1:s(2)
    for j=1:s(1)
        if nd == 1
            out(j, i) = in(j, i)-mean(in(:, i));
        else
            l = 1;
            for k = 1:length(squeeze(in(j, i, :)))
                out(j, i, l) = in(j, i, l) - mean(in(:, i, l));
                l=l+1;
            end
        end
    end
end