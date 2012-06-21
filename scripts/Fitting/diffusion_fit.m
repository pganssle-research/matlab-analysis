function f = diffusion_fit(x, G, n, tau)
    % Fits the data to a diffusion measurement.
    % Eqn: x(2)*exp(-((hg*tau*t).^2)*n*tau*x(1)/3);
    
    hg = 4257.77481*2*pi;
    
	 tau = tau*1e-3;
	 t = 2*n*tau;
	 D = x(1);
	 A = x(2);
	 
    f = A*exp(-((hg*tau*G).^2)*t*D/3);
   