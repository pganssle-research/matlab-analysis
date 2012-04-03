function [i_out] = dir2indir(curves, i_t)
    % Function which creates a matrix of the indirect dimensions from the
    % indirect dimensions.
    
    % Find out what we're working with
    % s(1) = # of direct points
    % s(2) = # of indirect points
    s = size(curves);
    
    