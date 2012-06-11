% Demonstration script for iPeak function. It generates a test signal
% consisting of several peaks, adds random noise, runs ipeak, then prints out 
% the actual and the measured peak positions, heights, widths and areas. 
% Each time you run it you get a different sample of noise. You can easily
% evaluate the accuracy of the measurements because the actual peak
% parameter values in this simulation are always integers.
%   T. C. O'Haver, September 2011
increment=1;
x=[1:increment:500];

% For each simulated peak, compute the amplitude, position, and width
amp=[1 2 3 4];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200 300 400];   % Positions of the peaks (Change if desired)
wid=[30 30 30 30];   % Widths of the peaks (Change if desired)
Noise=.5; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=30;  % Larger-than-usual values of SmoothWidth and FitWidth 
FitWidth=30; % will help deal with the very poor signal-to-noise ratio.

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('Detection and measurement of 4 peaks in a very noisy signal')
disp(['Signal-to-noise ratio of first peak = ' num2str(amp(1)/Noise) ] )
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks

if size(ActualPeaks)==size(MeasuredPeaks),
    PercentErrors=100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks
    AveragePercentErrors=mean(abs(100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks))
end
disp('The peak at x=100 is usually detected, but the accuracy of peak ')
disp('parameter measurement is poor because of the low signal-to-noise ratio.')
disp('Jump to the next/previous peaks using the Spacebar/Tab keys.')
