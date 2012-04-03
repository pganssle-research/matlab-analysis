function [curr, pow, curr2, pow2, powc] = calc_curr(R, R2, Rc)
    % Calculates the current across a coil (Rc) and the power across the
    % limiting resistor (R, pow), as well as the total current drawn
    % (curr2) and the power across the main limiting resistor (R2, pow2)
    %
    %                                    R -I->
    %                           +---/\/\/\/\/\---+
    %                     R2    |                |
    %        +5V ----/\/\/\/\--+              {Coil} (Rc) | (Pc)
    %                   -I2->   |        R       |         V
    %                           +---/\/\/\/\/\---+ 
    %
    %
    % [I, P, I2, P2, Pc] = calc_curr(R, R2, Rc);
    %  (Note, Ic = I)

    rtc = ((2*R+Rc)/(R^2+R*Rc))^(-1);  % Effective r.
    v1 = 5*(rtc/(rtc+R2));
    curr = v1/(Rc+R);
    pow = curr^2*R;
    
    curr2 = (5-v1)/R2;
    pow2 = curr2^2*R2;
    powc = curr^2*Rc;