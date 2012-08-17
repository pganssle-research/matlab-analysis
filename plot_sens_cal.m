pb = 'C:\\Omega\\Data\\Calibration\\CoilTest\\120720-Calibration-%03dHz-0000.mcd';

f = [2:4, 5:5:300, 310:10:490];

ftd = zeros(1024, length(f));
fs = zeros(1024, length(f));
as = zeros(1, length(f));


for i = 1:length(f)
	s = sprintf(pb, f(i));
	
	if(exist(s, 'file'))
		out = mc_read_data(s);
		
		fs(:, i) = out.f;
		ftd(:, i) = mean(out.fft, 2);
		
		hold all
		%plot_mc_struct(out, [-1], 2);
		hold off;
		
		[~, mini] = min(abs(out.f-(f(i)-10)));
		[~, maxi] = min(abs(out.f-(f(i)+10)));
		
		as(i) = sum(abs(out.fft(mini:maxi)));
		
	end
end

plot(f, as);

%plot(fs, abs(ftd));

