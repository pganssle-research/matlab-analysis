function f = sin_fit(x, t)
    % Fits the data to two exponential decays.
    % Takes 4 arguments
    %
    % Equation: x(1)*sin(x(2)*t + x(3)/pi()) + x(4)
    
    % It's just a decaying sine wave.
    f = x(1)*sin(x(2)*t + x(3)*pi()) + x(4);