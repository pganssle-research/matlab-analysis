function f = exponential_fit(x, t)
    % Fits the data to a single exponential decay.
    % Takes 3 arguments - x(2)*exp(-t/x(1)) + x(3)
    
    % It's just an exponential decay with some offset.
    f = x(2)*exp(-t/x(1));
   