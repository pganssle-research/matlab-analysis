function f = bi_exp_fit(x, t)
    % Fits the data to two exponential decays.
    % Takes 5 arguments - x(2)*exp(-t/x(1)) + x(4)*exp(-t/x(3))
    
    % It's just an exponential decay with some offset.
    f = x(2)*exp(-t/x(1)) + x(4)*exp(-t/x(3));
   