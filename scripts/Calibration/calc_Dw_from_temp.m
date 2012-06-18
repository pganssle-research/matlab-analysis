function Dw = calc_Dw_from_temp(temp)
    % Based on a calibration from Wang, J.H. J. Am. Chem. Soc. 1951, 73(2)
    % Plug in the temperature and it will calculate the self-diffusion
    % constant of water. Temp is in Celcius
    %
    % Answer is in cm^2/s.
    %
    % Dw = calc_Dw_from_temp(temp)
    
    temp = temp+273.15;
    
    % Viscosity will be approx A*10^(B/(T-C));
    A = 2.414e-4; % In poise
    B = 247.8; % In K
    C = 140; % In K
           
    Dws = [1.05, 1.26, 1.65, 1.77, 2.13, 2.5, 2.75, 3.48, 4.35]*1e-5;
    temps = [0, 5, 10, 17, 25, 28, 35, 45, 55]+273.15;
   
    visc = A*10.^(B./(temps-C));
  
    p = polyfit(temps, (Dws.*visc), 1);
    
    Dw = polyval(p, temp)./(A*10.^(B./(temp-C)));
   
end