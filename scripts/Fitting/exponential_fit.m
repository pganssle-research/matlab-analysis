function f = exponential_fit(x, t)
    % Fits the data to a sum of exponential decays. length(x) must be even.
    % A = x(2:2:end)
	 % tau = (x(1:2:end)
	 % f = sum(A.*exp(-t./tau));
    	 
    % It's just an exponential decay with some offset.
	 tau = x(1:2:end);
	 A = repmat(x(2:2:end), length(t), 1);
	 
	 f = sum(A.*exp(-t*(1./tau)), 2)';
	  
	 a = 1;