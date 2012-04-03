function f = exp_sin_fit(x, t)
    % Fits the data to two exponential decays.
    % Takes 5 arguments
    %
    % Equation: x(1)*exp(-t/x(2))*cos(x(3)*t + x(4)*pi) + x(5)
    
    % It's just a decaying sine wave.
    f = x(1)*exp(-t*2/x(2)).*cos(x(3)*t + x(4)*pi) + x(5);
   