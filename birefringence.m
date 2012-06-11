function out = birefringence(dphi, theta)
phi = [0; dphi];

phix = phi(1);
phiy = phi(2);

M = [exp(1i*phix), 0; 0, exp(1i*phiy)];
M = rot(theta)*(M*rot(-theta));

M = [exp(1i*phix)*(cos(theta).^2) + exp(1i*phiy)*(sin(theta).^2), (exp(1i*phix) - exp(1i*phiy)) * cos(theta)*sin(theta);
    (exp(1i*phix) - exp(1i*phiy)) * cos(theta)*sin(theta), exp(1i*phix)*(sin(theta).^2) + exp(1i*phiy)*(cos(theta).^2)];

out = M;

function out = rot(theta)
out = [cos(theta), -sin(theta); sin(theta), cos(theta)];
