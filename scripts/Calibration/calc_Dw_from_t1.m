function Dw = calc_Dw_from_t1(t1)
    % Uses calc_temp_from_T1 and calc_Dw_from_temp to calculate the
    % diffusion constant based on the T1.
    
    Dw = calc_Dw_from_temp(calc_temp_from_T1(t1));
   
end