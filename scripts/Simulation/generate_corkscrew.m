% Generate a corkscrew
cork = zeros(100, 100);
cork2 = zeros(100, 100);
for i = 1:100
    roll = zeros(10, 10);
    roll(round(4*cos((i/10)*pi))+6, round(4*sin((i/10)*pi))+6) = 1;
    
    cork(:, i) = reshape(roll, 100, 1);
    cork2(:, i) = reshape(roll', 100, 1);
end