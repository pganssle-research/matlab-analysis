function out = simulate_pulse_error(e, field, gap, np)
% Simulates the output of pulses with error.
% e is the error vector in percent ([x, y, z] or scalar)
% field is the field, in microgauss (default = 0)
% gap is the interpulse spacing in ms (default: 60)
% np is number of pulses (default: 100)
%
% Output sampling rate is 1ks/S;
%
% out = simulate_pulse_error(e[, field, gap, np);

if(~exist('e', 'var') || (~isscalar(e) && length(e) ~= 3))
   error('Need an error vector'); 
end

if(isscalar(e))
    e = repmat(e, 3, 1);
end

if(~exist('field', 'var'))
   field = [0, 0, 0]; 
end

if(~exist('gap', 'var'))
   gap = 60; 
end

if(~exist('np', 'var'))
   np = 100; 
end

out = zeros(3, np*gap);
x = [1; 0; 0];
y = [0; 1; 0];
z = [0; 0; 1];

pulses = RotMat(x, pi/2-((pi/2)*e(1)/100))*RotMat(y, pi-pi*e(2)/100)*RotMat(x, pi/2-(pi/2*e(1)/100));
%pulses = RotMat(x, pi-pi*e(1)/100);

gH = 4.257*1e-6; % Gyromagnetic ratio in kHz/uG

fangle = 2*pi*norm(field)*gH; % Angle in, I believe, radians.
f = RotMat(field, fangle);

vec = z;

for i = 1:np
    j = (i-1)*gap;
    vec = pulses*vec;
    
    if(fangle == 0)
        out(:, (j+1):(j+gap)) = repmat(vec(:), 1, gap);
    else
        for k = 1:gap
           vec = f*vec;
           out(:, j+k) = vec(:);
        end
    end
end

function mat = RotMat(axis, angle)
if(size(axis, 2) > size(axis, 1))
    axis = axis';
end
I = diag(ones(3, 1));
axis = axis/norm(axis);
cp = cos(angle);
ax = [0, -axis(3), axis(2); axis(3), 0, -axis(1); -axis(2), axis(1), 0];

mat = I*cp + (1-cp)*(axis*axis') + ax*sin(angle);