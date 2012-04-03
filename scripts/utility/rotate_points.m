function points = rotate_points(p, a1, a2, thet)
% Rotates points about axis defined by a1 -> a2 angle thet.
l1 = length(a1);
l2 = length(a2);

if(l1 ~= l2)
    error('Axis points must be the same size.');
elseif(l1 > 3)
    error('Rotations in higher dimensions not yet supported.')
else
    if(size(a1, 1) ~= 1)
        a1 = a1';
    end
    
    if(size(a2, 1) ~= 1)
        a2 = a2';
    end
    
    a1 = zeros(1, 3) + a1;
    a2 = zeros(1, 3) + a2;
end

t = false;
if(size(p, 2) > 3)
    p = p';
    t = true;
end

lp = size(p);

if(size(p, 2) == 1)
    p(:, 2:3) = 0;
elseif(size(p, 2) == 2)
    p(:, 3) = 0;
end

for i = 1:3
    p(:, i) = p(:, i)-a1(i);
end

R = R3d(thet, a2-a1);

for i = 1:size(p, 1)
    p(i, :) = (R*p(i, :)')';
end

for i = 1:3
    p(:, i) = p(:, i)+a1(i);
end

if(lp(2) == 2)
    p(:, 3) = [];
elseif(lp(2)==1)
    p(:, 2:3) = [];
end
    
if(t)
    p = p';
end

points = p;

function R=R3d(deg,u)
%R3D - 3D Rotation matrix counter-clockwise about an axis.
%
%R=R3d(deg,axis)
%
%Input is in degrees.
%
%See also Rx,Ry,Rz,R3d,M2d,M3d

R=eye(3);
u=u(:)/norm(u);
x=deg; %abbreviation

for ii=1:3
    v=R(:,ii);
    
    R(:,ii)=v*cosd(x) + cross(u,v)*sind(x) + (u.'*v)*(1-cosd(x))*u;
      %Rodrigues' formula
end