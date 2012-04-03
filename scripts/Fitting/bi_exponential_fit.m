function f = bi_exponential_fit(x, t)
    % Fits the data to a single exponential decay.
    % Takes 5 arguments - x(2)*exp(-t/x(1)) + x(3)*exp(-t/x(4)) + x(5)
    
    % It's just an exponential decay with some offset.
    f = x(2)*exp(-t/x(1)) + x(3);
   