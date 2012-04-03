function T = calc_temp_from_T1(t1)
    % Based on the calibration done by Scott on the Magritek, this
    % calculates the temperature of the sample based on the t1.
    %
    % T = calc_temp_from_T1(t1)
    
    Ts = [27,  39, 54, 64, 84];
    T1s = [2.5, 3.2, 4.5, 5.7, 7.8];
    
    p = polyfit(T1s, Ts, 1);
    
    T = polyval(p, t1);
end