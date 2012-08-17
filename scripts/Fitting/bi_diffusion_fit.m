function f = bi_diffusion_fit(x, G, n, tau)
    % Fits the data to a diffusion measurement.
    % Eqn: x(2)*exp(-((hg*tau*t).^2)*n*tau*x(1)/3) 
	 %    + x(4)*exp(-((hg*tau*t).^2)*n*tau*x(3)/3);
    
    hg = 4257.77481*2*pi;
    
	 tau = tau*1e-3;
	 t = 2*n*tau;
	 D1 = x(1);
	 A1 = x(2);
	 D2 = x(3);
	 A2 = x(4);
	 
    f = A1*exp(-((hg*tau*G).^2)*t*D1/3) + A2*exp(-((hg*tau*G).^2)*t*D2/3);
   