function T = calc_temp_from_T1(t1, t1_or_t2)
    % Based on the calibration done by Scott on the Magritek, this
    % calculates the temperature of the sample based on the t1.
    %
    % T = calc_temp_from_T1(t1)
	 
	 if(~exist('t1_or_t2', 'var'))
		 t1_or_t2 = 1;
	 end
	 
	 if(t1_or_t2 == 2)
		 T1s = [2.24, 3.19, 4.9, 6.2, 6.7];
	 else
		T1s = [2.5, 3.2, 4.5, 5.7, 7.8];
	 end
	 
    Ts = [27, 39, 54, 64, 84];
    
    p = polyfit(T1s, Ts, 1);
    
    T = polyval(p, t1);
end