function out = get_peaks(freq, mag)
	% Function for getting peaks
	% 
	% Usage: out = get_peaks(freq, mag)
	%
	%
	% Must be only 2D
	% Grabs the top 3 peaks in each spectrum and pumps them into a
	% number_of_spectra*2*3 array
	% out(i, 1, 1) is the frequency of the biggest peak in the ith spectrum
    % out(i, 2, 1) is the magnitude of that peak
    % out(i, 1, 2) is the frequency of the 2nd biggest peak, etc.	

	s = size(mag);
	out = zeros(s(2), 2, 3);
	for i = 1:s(2)
		p = ipeak([freq mag(:, i)]);
		[z, ix] = sort(p(:, 3), 'descend');
		
		if i == 1 || i == 5
			z(1:3)
		end
		
		out(i, 1, 1:3) = p(ix(1:3), 2);
		out(i, 2, 1:3) = z(1:3);
	end