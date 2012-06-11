% Demonstration script for iPeak function. It generates a test signal
% consisting of several peaks, adds random noise, runs ipeak, then prints out 
% the actual and the measured peak positions, heights, widths and areas. 
% Each time you run it you get a different set of peaks. You can easily
% evaluate the accuracy of the measurements because the actual peak
% parameter values in this simulation are always integers.
%   T. C. O'Haver, September 2011
increment=1;
x=[1:increment:4000];

% For each simulated peak, compute the amplitude, position, and width
amp=round(5.*randn(1,38));  % Amplitudes of the peaks  (Change if desired)
pos=[200:100:3900];   % Positions of the peaks (Change if desired)
wid=30.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=0.1; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.2,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=round(WidthPoints/2);  % SmoothWidth should be roughly equal to 1/2 the peak width (in points)
FitWidth=round(WidthPoints/2); % FitWidth should be roughly equal to 1/2 the peak widths(in points)

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW,
% with AUTOZERO=1 (ON) and the top window showing the first peak. Once
% iPeak is running, you can change any of these using 
% keystroke commands (press "K" to see a list of commands).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),120,1);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks


if size(ActualPeaks)==size(MeasuredPeaks),
    PercentErrors=100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks
    AveragePercentErrors=mean(abs(100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks))
end

disp('Demonstration of autozero background correction, for separated, narrow ')
disp('peaks on a large baseline. Hint: Turn on the Autozero mode (T key),')
disp('adjust the zoom setting so that the peaks are shown one at a time in')
disp('the upper window, then press the P key to display the peak table.')
disp('Jump to the next/previous peaks using the Spacebar/Tab keys.')

