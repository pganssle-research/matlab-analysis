function f = exponential_fit(x, t)
    % Fits the data to a sigmoid function.
    % Takes 6 arguments: f = x(1) ./ (x(2) + x(3)*exp(-(x(4)*t - x(5)))) +
    % x(6)
    
    f = x(1) ./ (x(2) + x(3)*exp(-(x(4)*t - x(5)))) + x(6);
   