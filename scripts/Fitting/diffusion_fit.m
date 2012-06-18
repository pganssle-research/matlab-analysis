function f = diffusion_fit(x, t)
    % Fits the data to a diffusion measurement.
    % Eqn: x(2)*exp(-((hg*tau*t).^2)*n*tau*x(1)/3);
    
    hg = 4257*2*pi;
    
    n = 24;
    tau = 60e-3;
    
    f = x(2)*exp(-((hg*tau*t).^2)*n*tau*x(1)/3);
   