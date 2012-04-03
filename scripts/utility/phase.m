function out = phase(data, thet, rad)
% Gives a real phased spectrum, thet is normally in degrees, but passing 1
% to rad converts from radians (i.e. thet can be given in radians).
%
% out = phase(data, thet, rads);

if(nargin < 3)
    rad = 0;
end

% If it's in radians, convert to degrees
if(rad)
    thet = thet*180/pi;
end

st = size(thet);
if(st(2) > 1)
    thet = thet';
end
st = size(thet);

thet = thet*ones(size(data))';
data = ones(st)*data';

out = cos(thet).*real(data) +sin(thet).*imag(data);
end