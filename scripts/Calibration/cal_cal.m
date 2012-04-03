function [field, response] = cal_calc(pulse_length, voltage, resistance)
	% Calculate the response of coils based on 180 pulses on hydrogen spins.
	% pulse = time in microseconds
	% voltage = peak voltage in volts
	% resistance = resistance in ohms
	% 
	% When a pulse is 180, then t = (2n+1)*pi/omega
	% So assuming n = 0, omega = pi/t.
	% Convert to frequency units, f = 1/2t
	% hgam*B = f, B = 1/(2*hgam*t)
	% Since we produced the field with a pulse, the response is just B/I.

	amps = voltage/resistance; % Current in amps
	p_t = pulse_length/1e6; % Pulse time in seconds
	hgam = 4258; % Gyromagnetic ratio in Hz/G
	
	field = 1/(2*hgam*p_t); % Absolute field in Gauss
	response = field/amps; % Response in Gauss/Amp