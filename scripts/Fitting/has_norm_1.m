function out = has_norm_1(y)
arrayfun(@(x)norm(y(x:(x+2))), 1:3:length(y))
