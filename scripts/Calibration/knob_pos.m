function pos = knob_pos(voltage, x_y_z)
	% What should the knob position be to get a certain voltage?
    %
    % Usage: pos = knob_pos(voltage, x_y_z)
    %
    % x_y_z = 1 for x, 2 for y, 3 for z
    % Voltage should be in millivolts.

	if x_y_z == 1
        p = [2.6024 -14.1324];
    elseif x_y_z == 2
        p = [2.5205 -13.7850];
    elseif x_y_z == 3
        p = [1.4197 -8.3862];
    else
        pos = -1;
        return
    end
    
    pos = (voltage - p(2))/p(1);