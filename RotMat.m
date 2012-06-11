function mat = RotMat(axis, angle)
if(size(axis, 2) > size(axis, 1))
    axis = axis';
end
I = diag(ones(3, 1));
axis = axis/norm(axis);
cp = cos(angle);
ax = [0, -axis(3), axis(2); axis(3), 0, -axis(1); -axis(2), axis(1), 0];

mat = I*cp + (1-cp)*(axis*axis') + ax*sin(angle);