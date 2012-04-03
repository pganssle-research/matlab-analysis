data = [mice(:).i_ref];
uc = unique(cols);
od = cell(length(uc), 1);

figure;
hold on
for i = 1:length(uc)
    temp = logical(cols == uc(i));
    temp = ones(sum(temp), 2);
    temp(:, 1) = data(logical(cols == uc(i)));
    temp(:, 2) = temp(:, 2)*uc(i);
    od{i} = temp;
    
    h = scatter(temp(:, 1), temp(:, 2), 30, temp(:, 2));
end 
hold off