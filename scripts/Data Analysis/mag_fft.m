function [out, f] = mag_fft(data, sr)
	% Takes data, returns magnitude fourier transform and the frequency vector.
	%
	% Usage: [out, f] = mag_fft(data);

	s = size(data);
	
	% Use the next highest power of 2 greater than or equal to the length to calculate the fft
	nfft = 2^(nextpow2(s(1)));
	
	% Calculate number of unique points.
	n_u_points = ceil((nfft+1)/2);
	
	if length(s) > 1
		i = 1;
		out = zeros([n_u_points s(2:end)]);
		for d = data
			% This is a weird way to do this, so I'll explain.
			% Out is a matrix of indeterminite dimension, so we're using a "foreach" loop
			% which will give us the column in the first dimension for each subsequent
			% dimension. This of course doesn't give us the right index
			out(:, i) = mag_fft_1d(d);
			i = i+1;
		end
		
	else
		out = mag_fft_1d(sr);
	end
	
	f = (0:n_u_points-1)*sr/nfft;

end
function [s] = mag_fft_1d(d)
	% Subfunction for doing the ffts 1 at a time.
	% Only feed this vectors.
	% This is basically directly ripped off from here: http://www.mathworks.com/support/tech-notes/1700/1702.html
	
	% Use the next highest power of 2 greater than or equal to the length to calculate the fft
	nfft = 2^(nextpow2(length(d)));
	
	% Calculate number of unique points.
	n_u_points = ceil((nfft+1)/2);
	
	% FFT with zero-padding.
	fftd = fft(d);
	
	% FFT is symmetric, throw away the second half.
	fftd = fftd(1:n_u_points);
	
	% Magnitude, scaled appropriately, not a function of length of x.
	m_fftd = abs(fftd)/length(d);
	%m_fftd = sqrt(real(fftd).^2 + imag(fftd).^2);
	
	% Square it.
	 m_fftd = m_fftd.^2;
	
	% We dropped half the FFT, so we multiply it by 2 to keep the same energy.
	% The DC component and the nyquist component are unique and should not be multiplied by 2.
	
	if rem(nfft, 2) % odd nfft exclude Nyquist point.
		m_fftd(2:end) = m_fftd(2:end)*2;
	else
		m_fftd(2:end-1) = m_fftd(2:end-1)*2;
	end
	
	s = m_fftd;
end