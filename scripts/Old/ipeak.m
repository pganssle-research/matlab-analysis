function P=ipeak(DataMatrix,PeakD,AmpT,SlopeT,SmoothW,FitW)
% ipeak(DataMatrix)
% Keyboard-operated Interactive Peak Finder for data in a single data
%  matrix "DataMatrix", with x values in row 1 and y values in row 2.
% Returns the peak table in P (Peak #, Position, Height, Width)
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% T. C. O'Haver (toh@umd.edu), Version 1.2, April 2011.  
%  
% EXAMPLE 1:   x=[0:.1:100];
%              y=(x.*sin(x)).^2;
%              ipeak([x' y']);
% Displays DataMatrix [x' y'] with arbitrary starting
% values for peak parameters AmpT, SlopeT, SmoothW, FitW.
%
% EXAMPLE 2: x=[0:.1:100];
%            y=5+5.*cos(x)+randn(size(x));
%            ipeak([x' y'],10);
% The additional argument - 10 in this case - is an estimate of the peak
% density (maximum number of peaks that would fit into the data record).
%
% EXAMPLE 3: ipeak([x' y'],0,.0001,33,33);
% As above, but specifies initial values of AmpT, SlopeT, SmoothW, FitW.
%
% Keyboard Controls:
% Pan signal left and right:  Coarse pan: < and >   
%                             Fine pan: left and right cursor arrow keys
% Zoom in and out:            Coarse zoom: / and '   
%                             Fine zoom: up and down cursor arrow keys
% Adjust AmpThreshold:        A,Z  (A increases, Z decreases)
% Adjust SlopeThreshold:      S,X  (S increases, X decreases)
% Adjust SmoothWidth:         D,C  (D increases, C decreases)
% Adjust FitWidth:            F,V  (F increases, Z decreases)
% Baseline:                   B, then click baseline at 8 points
% Print peak table:           P  Prints Peak #, Position, Height, Width
% Print keyboard commands:    K  Prints this list
% Print parameters            Q (Prints AmpT, SlopeT, SmoothW, FitW)
% Print report                R  Prints Peak table and parameters
% Peak labels ON/OFF:         L
%
% For large data sets, to view only a portion of the data over a 
% specified x-axis range, you can type 
% n1=val2ind(x,x1);n2=val2ind(x,x2);ipeak([x(n1:n2)' y(n1:n2)'])
% where x1 and x2 are the end limits of the x-axis values displayed.
global X
global Y
global xo
global dx
global SlopeThreshold 
global AmpThreshold  
global SmoothWidth
global FitWidth
global PeakLabels
global P
close
format short g
format compact
warning off all
% If DataMatrix is in the wrong transposition, fix it.
datasize=size(DataMatrix);
if datasize(1)<datasize(2),DataMatrix=DataMatrix';end
% Assign arguments to internal global variables
X=DataMatrix(:,1); % Split matrix argument 
Y=DataMatrix(:,2);
switch nargin
    % 'nargin' is the number of arguments
    case 1
    % Calculate default values of peak detection parameters
      PeakDensity=50;   
      % Estimate approx number of points in a peak half-width
      WidthPoints=length(Y)/PeakDensity;  
      SlopeThreshold=WidthPoints^-2;  
      AmpThreshold=abs(min(Y)+0.02*(max(Y)-min(Y))); 
      SmoothWidth=round(WidthPoints/3);  
      FitWidth=round(WidthPoints/3); 
  case 2            
     % Calculate values of peak detection parameters
     % arguments based on the peak density, PeakD
      PeakDensity=PeakD;    
      % Estimate approx number of points in a peak half-width
      WidthPoints=length(Y)/PeakDensity;  
      SlopeThreshold=WidthPoints^-2;  
      AmpThreshold=abs(min(Y)+0.02*(max(Y)-min(Y)));  
      SmoothWidth=round(WidthPoints/3);  
      FitWidth=round(WidthPoints/3); 
  case 6   
      % Initial values of all peak detection parameters specified in arguments 
      PeakDensity=PeakD;  
      SlopeThreshold=SlopeT;
      AmpThreshold=AmpT;
      SmoothWidth=SmoothW;
      FitWidth=FitW; 
  otherwise
      disp('Invalid number of arguments')
end % switch nargin

if FitWidth<3,FitWidth=3;end   % Keep FitWidth above 2 
xo=length(Y)/2; % Initial Pan setting
dx=length(Y)/4; % Initial Zoom setting
PeakLabels=0; % Start with onlt  numbering peaks

% Plot the signal 
%RedrawSignal(X,Y,xo,dx);

% Attaches KeyPress test function to the figure.
%set(gcf,'KeyPressFcn',@ReadKey)
%uicontrol('Style','text')
% end of outer function
% ----------------------------SUBFUNCTIONS--------------------------------
function ReadKey(obj,eventdata)
% Interprets key presses from the Figure window.
% When a key is pressed, interprets key and calls corresponding function.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
global X
global Y
global xx
global yy
global xo
global dx
global SlopeThreshold 
global AmpThreshold  
global SmoothWidth
global FitWidth
global PeakLabels
global P
key=get(gcf,'CurrentCharacter');
if ischar(key),
  switch double(key),
    case 29
        % Pans one point down when left arrow pressed.
          xo=xo-1;
          if xo<1,xo=1;,end
          [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 28
        % Pans one point up when right arrow pressed.
        ly=length(Y);
        xo=xo+1;
        if xo>ly,xo=ly;,end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 46
        % Pans 2% down when < key pressed.
        ly=length(Y);
        xo=xo-ly/50;
        if xo<1,xo=1;,end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 44
        % Pans 2% up when > key pressed.
        ly=length(Y);
        xo=xo+ly/50;
        if xo>ly,xo=ly;,end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 31
        % Zooms one point up when up arrow pressed.
        dx=dx+2;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 30
        % Zooms one point down when down arrow pressed.
        dx=dx-2;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 47
        % Zooms 2% up when / pressed.
        ly=length(Y);
        dx=dx+ly/50;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);
    case 39
        % Zooms 2% down when ' pressed.
        ly=length(Y);
        dx=dx-ly/50;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);      
    case 98
        % When 'b' key is pressed, user clicks graph 
        % to enter background points, then graph re-drawn.
        BaselinePoints=8;  % Change as you wish <============
        % Acquire background points from user mouse clicks
        subplot(2,1,2)
        title(['Click on ' num2str(BaselinePoints) ' points on the baseline between the peaks.'])
        bX=[];bY=[];
        for g=1:BaselinePoints;
           [clickX,clickY] = ginput(1);
           bX(g)=clickX;
           bY(g)=clickY;
           xlabel(['Baseline point '  num2str(g) ' / ' num2str(BaselinePoints) ])
        end
        yy=Y;
        for k=1:length(bX)-1,
           fp=val2ind(X,bX(k)); % First point in segment
           lp=val2ind(X,bX(k+1));  % Last point in segment
           % Subtract piecewise linear background from Y
           yy(fp:lp)=Y(fp:lp)-((bY(k+1)-bY(k))/(bX(k+1)-bX(k))*(X(fp:lp)-bX(k))+bY(k));
        end
        Y=yy;
        % Estimate initial value of AmpThreshold
        AmpThreshold=min(Y)+0.02*(max(Y)-min(Y));  
        [xx,yy]=RedrawSignal(X,Y,xo,dx);    
    case 97
        % When 'a' key is pressed, increases "AmpThreshold" by 10%
        AmpThreshold=abs(AmpThreshold+.1*AmpThreshold);
        [xx,yy]=RedrawSignal(X,Y,xo,dx);      
    case 122
        % When 'z' key is pressed, decreases "AmpThreshold" by 10%
        AmpThreshold=AmpThreshold-.1*AmpThreshold;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);      
   case 115 % When 's' key is pressed, increases "SlopeThreshold" by 10%
         SlopeThreshold=SlopeThreshold+.1*SlopeThreshold;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
   case 120 % When 'x' key is pressed, decreases "SlopeThreshold" by 10%
         SlopeThreshold=SlopeThreshold-.1*SlopeThreshold;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
   case 100
        % When 'd' key is pressed, increases "SmoothWidth" by 1
        SmoothWidth=SmoothWidth+1;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 99
        % When 'c' key is pressed, decreases "SmoothWidth" by 1
        SmoothWidth=SmoothWidth-1;
        if SmoothWidth<1, SmoothWidth=1;,end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
     case 102
        % When 'f' key is pressed, increases "FitWidth" by 1
        FitWidth=FitWidth+1;
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 118
        % When 'v' key is pressed, decreases "FitWidth" by 1
         FitWidth=FitWidth-1;
         if FitWidth<3, FitWidth=3;,end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);   
    case 114
        % When 'r' key is pressed, prints
        disp('--------------------------------------------------------')
        disp(['Amplitude Threshold (AmpT) = ' num2str(AmpThreshold) ] )
        disp(['Slope Threshold (SlopeT) = ' num2str(SlopeThreshold) ] )
        disp(['Smooth Width (SmoothW) = ' num2str(SmoothWidth) ] )
        disp(['Fit Width (FitW) = ' num2str(FitWidth) ] )
        disp('         Peak #    Position      Height       Width')
        disp(P)
    case 112
        % When 'p' key is pressed, prints out peak table 
        disp('         Peak #    Position      Height       Width')
        disp(P)
    case 107
        % When 'k' key is pressed, prints out table of keyboard commands
        disp('Keyboard Controls:')
        disp(' Pan signal left and right:  Coarse pan: < and >')   
        disp('                             Fine pan: left and right cursor arrows')
        disp(' Zoom in and out:            Coarse zoom: / and "  ') 
        disp('                             Fine zoom: up and down cursor arrows')
        disp(' Adjust AmpThreshold:        A,Z')
        disp(' Adjust SlopeThreshold:      S,X')
        disp(' Adjust SmoothWidth:         D,C')
        disp(' Adjust FitWidth:            F,V')
        disp(' Baseline:                   B, then click baseline at 8 points')
        disp(' Print report                R')
        disp(' Print peak table:           P')
        disp(' Print keyboard commands:    K')
        disp(' Print parameters:           Q' )
        disp(' Peak labels ON/OFF:         L' )
     case 113
        % When 'Q' is pressed, prints peak detection parameters on a single line     
        disp([ num2str(AmpThreshold) ',' num2str(SlopeThreshold)  ',' num2str(SmoothWidth)  ',' num2str(FitWidth) ])
    case 108
        % When 'L' is pressed, toggles on/off peak labels in upper panel
        if PeakLabels==0,
            PeakLabels=1;
        else
            PeakLabels=0;
        end
        [xx,yy]=RedrawSignal(X,Y,xo,dx);  
    otherwise  
       UnassignedKey=double(key)
       disp('Press k to print out list of keyboard commands')
   end % switch
end % if
% ----------------------------------------------------------------------    
function [xx,yy]=RedrawSignal(X,Y,xo,dx);
% Plots the entire signal (X,Y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys.
global SlopeThreshold 
global AmpThreshold  
global SmoothWidth
global FitWidth
global PeakLabels
global P
Startx=round(xo-(dx/2));
Endx=round(xo+(dx/2)-1);
if Endx>length(Y),
    Endx=length(Y);
end
if Startx<1,
     Startx=1;
end
PlotRange=[Startx:Endx];
if PlotRange<5, PlotRange=[xo:xo+5];end
xx=X(PlotRange);
yy=Y(PlotRange); 
hold off
% clf
% Plots isolated segment (xx,yy) in the upper half
figure(1);subplot(2,1,1);plot(xx,yy,'g.'); 
hold on
title('Pan and Zoom to inspect peaks. Press K for keyboard commands')
axis([X(Startx) X(Endx) min(yy) max(yy)]);

% Bottom half of the figure shows full signal
figure(1);subplot(2,1,2);cla
plot(X,Y,'g')  % Graph the signal
axis([X(1) X(length(X)) min(Y) max(Y)]); % Update plot
P=findpeaks(X,Y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
title(' Adjustment Keys   AmpT: A/Z    SlopeT: S/X   SmoothWidth: D/C    FitWidth: F/V')
xlabel(['    AmpT = ' num2str(AmpThreshold) '       SlopeT = ' num2str(SlopeThreshold) '    SmoothW = ' num2str(SmoothWidth) '    FitW = ' num2str(FitWidth) ])
text(P(:,2),P(:,3),num2str(P(:,1)))  % Number the peaks found on the lower graph
hold on
% Mark the zoom range on the full signal with two magenta dotted vertical lines
plot([min(xx) min(xx)],[min(Y) max(Y)],'m--')
plot([max(xx) max(xx)],[min(Y) max(Y)],'m--') 
% Number the peaks found on the upper graph
subplot(2,1,1);
if PeakLabels==1,
   % Label the peaks on the upper graph with number, position, height, and
   % width
   topaxis=axis;
   yrange=topaxis(4)-topaxis(3);
   pos1=.1*yrange;
   pos2=.2*yrange;
   pos3=.3*yrange;
   pos4=.4*yrange;
   PP=findpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
   text(P(:,2),P(:,3)-pos1,num2str(P(:,1)))
   text(PP(:,2),PP(:,3)-pos2,num2str(PP(:,2)))
   text(PP(:,2),PP(:,3)-pos3,num2str(PP(:,3)))
   text(PP(:,2),PP(:,3)-pos4,num2str(PP(:,4)))
else   
   topaxis=axis;
   yrange=topaxis(4)-topaxis(3);
   pos1=.1*yrange; 
   % Number the peaks on the upper graph
   text(P(:,2),P(:,3)-pos1,num2str(P(:,1))) 
end
markpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
hold off
% ----------------------------------------------------------------------
function markpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
% Draws blue best-fit line through data points near peak and computes
% peak position, amplitude, and width
PP=[0 0 0 0];
PP=findpeaks(xx,yy,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth);
hold on
if PP(1)>0,
   sizePP=size(PP);
   lengthPP=sizePP(1);
  for PeakNumber=1:lengthPP,
    FitWidth=round(real(FitWidth));
    subplot(2,1,1);
    if PeakNumber>lengthPP,PeakNumber=lengthPP,end
    n1=val2ind(xx,PP(PeakNumber,2))-round(FitWidth/2);
    n2=val2ind(xx,PP(PeakNumber,2))+round(FitWidth/2);
    if n1<1, n1=1;,end
    if n2>length(yy), n2=length(yy);,end
    PlotRange=[n1:n2];
    xxx=xx(PlotRange);
    yyy=yy(PlotRange);
    % Fit parabola to log10 of sub-group
    [coef,S,MU]=polyfit(xxx,log(abs(yyy)),2);  
    c1=coef(3);c2=coef(2);c3=coef(1);
    % Compute peak position and height or fitted parabola
    PeakX=-((MU(2).*c2/(2*c3))-MU(1));  
    PeakY=exp(c1-c3*(c2/(2*c3))^2);
    MeasuredWidth=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));
    subplot(2,1,1);
    plot(xxx,PeakY.*gaussian(xxx,PeakX,MeasuredWidth),'b');
  end  
end
% ----------------------------------------------------------------------
function P=findpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup)
% function P=findpeaks(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup)
% Function to locate the positive peaks in a noisy x-y data
% set.  Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number and
% position, height, and width of each peak. SlopeThreshold,
% AmpThreshold, and smoothwidth control sensitivity
% Higher values will neglect smaller features. Peakgroup
% is the number of points around the "top part" of the peak.
% T. C. O'Haver, 1995.  Version 2  Last revised Oct 27, 2006
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
d=fastsmooth(deriv(y),smoothwidth);
n=round(peakgroup/2+1);
P=[0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=smoothwidth:length(y)-smoothwidth,
   if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
       % if slope of derivative is larger than SlopeThreshold
     if d(j)-d(j+1) > SlopeThreshold*y(j), 
        if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
          for k=1:peakgroup, % Create sub-group of points near peak
              groupindex=j+k-n+1;
              if groupindex<1, groupindex=1;end
              if groupindex>vectorlength, groupindex=vectorlength;end
            xx(k)=x(groupindex);yy(k)=y(groupindex);
          end% Fit parabola to log10 of sub-group with centering and scaling
          [coef,S,MU]=polyfit(xx,log(abs(yy)),2);  
          c1=coef(3);c2=coef(2);c3=coef(1);
          % Compute peak position and height of fitted parabola
          PeakX=-((MU(2).*c2/(2*c3))-MU(1));   
          PeakY=exp(c1-c3*(c2/(2*c3))^2);
          MeasuredWidth=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));
          % if the peak is too narrow for least-squares technique to work
          % well, just use the max value of y in the sub-group of points near peak.
          if peakgroup<7,
             PeakY=max(yy);
             pindex=val2ind(yy,PeakY);
             PeakX=xx(pindex(1));
          end         
          % Construct matrix P. One row for each peak 
          % detected, containing the peak number, peak 
          % position (x-value) and peak height (y-value).
          P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth];
          peak=peak+1;
        end
      end
   end
end
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6006.*wid)) .^2);
% ----------------------------------------------------------------------
function SmoothY=fastsmooth(Y,smoothwidth)
%  fastsmooth(Y,w) smooths vector Y by triangular
% smooth of width = smoothwidth. Works well with signals up to 
% 100,000 points in length and smooth widths up to 1000 points. 
% Faster than tsmooth for smooth widths above 600 points.
% Example: fastsmooth([0 0 0 0 9 0 0 0 0],3) yields [0 0 1 2 3 2 1 0 0]
%  T. C. O'Haver, 2006.
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
for k=1:length(Y)-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s=s./w;
SumPoints=sum(s(1:w));
SmoothY=zeros(size(s));
for k=1:length(s)-w,
   SmoothY(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-s(k);
   SumPoints=SumPoints+s(k+w);
end
SmoothY=SmoothY./w;
% ----------------------------------------------------------------------
function d=deriv(a)
% First derivative of vector using 2-point central difference.
% Example: deriv([1 1 1 2 3 4]) yeilds [0 0 .5 1 1 1]
%  T. C. O'Haver, 1988.
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end
% ----------------------------------------------------------------------