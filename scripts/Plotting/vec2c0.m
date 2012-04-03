function out = vec2c0(vec)
% Function that gives you a cell of the form: {[0, vec(1)] [1, vec(2)], [2,
% 0]}
%
% For use with arrow3d, I guess.
% Usage: out = vec2c0(vec)

v1 = [0, vec(1)];
v2 = [0, vec(2)];
v3 = [0, vec(3)];

out = {v1, v2, v3};