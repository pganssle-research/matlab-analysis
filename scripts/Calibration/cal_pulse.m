function time = cal_resp(resp, voltage, resistance)
	% Returns the pi pulse length in microseconds
	%
	% resp = response in Gauss/Amp
	% voltage = peak voltage in volts
	% resistance = resistance in ohms
	% 
	% When a pulse is 180, then t = (2n+1)*pi/omega
	% So assuming n = 0, omega = pi/t.
	% Convert to frequency units, f = 1/2t
	% 
	% hgam*B = f (frequency units)
	% f = 1/(2t)
	% t = 1/(2*hgam*B)
    %
    % time = cal_resp(resp, voltage, resistance);
    % time = cal_resp(resp, current);

    if(nargin < 3)
        amps = voltage;
    else
        amps = voltage/resistance; % Current in amps
    end
    
    field = amps*resp; % Field
	hgam = 4258; % Gyromagnetic ratio in Hz/G
	
	time = 1/(2*hgam*field); % Response time in seconds
	time = time*1e6; % Response time in microseconds.
	
	