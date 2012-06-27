function f = diffusion_fit(x, G, n, tau)
    % Fits the data to a sum of diffusion fits.
	 %
	 % Amplitudes are the even components of (x)
	 % Diffusion coefficients are the odd components of (x)
	 % 
	 % f = sum(A.*exp(-((hg*tau*G)).^2*t*D/3));
	 
    hg = 4257.77481*2*pi;
    tau = tau*1e-3;
	 t = 2*n*tau;
	 
	 D = x(1:2:end);
	 A = repmat(x(2:2:end), length(G), 1);
	 
    f = sum(A.*exp(-((hg*tau*G').^2)*t*D/3), 2)';
   