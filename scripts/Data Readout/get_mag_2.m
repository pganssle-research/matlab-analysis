function [out, spans, points, times] = get_mag_2(start, len, num_windows, sr, window, data, asym, chan)
% Function for getting the magnetization curve as a function of time from an
% experiment consisting of a train of 180s.
%
% start = Time point to start from (ms)
% len = length of the plateaus (ms)
% num_windows = number of 'points' to acquire
% window = how long to average between pulses (ms)
% sr is the sampling rate
% data = the experimental data -> Always put in data with a vector size of [np, nchans, {maxsteps}]]
% Alternatively, pass this a struct and it will find the appropriate data.
% asym is an optional argument that allows for asymmetrical sampling of the medians
% chan is the channel you'd like to select (default is 1. 
%
% out is an output matrix of size (number of indirect points)*finish
%
%   Illustration of "start", "len" and "window", all in ms. Sample is taken
%   in the |__| regions.
%
%    ______        _|__|_        _|__|_
%   |      |      |      |      |      |
%   |      |      |      |      |      |
%   |      |_|__|_|      |_|__|_|      |
%   start -->|    |< len>|  <>
%                            ^-window
%
% Usage:
% [out, spans, points, times] = get_mag(start, freq, num_windows, sr, window, data, asym);

% Defaults
if nargin < 6
    fprintf('Insufficient number of arguments\n')
    return
end

if nargin < 7
    asym = 0.5;
end

if nargin < 8
    chan = 1;        
end

if(isstruct(data))
   if(isfield(data, 'adata'))
       data = data.adata;
   else
       data = data.mdata;
   end
end

% Select the channel
s = size(data);
if(length(s) < 2)
    return;
end

if(s(2) < chan)
    chan = 1;
    printf('Invalid channel selection, defaulting to 1.');
end

% Generate an equation to select the channel. Don't know a better way to do this. 
eqn = sprintf('squeeze(x(:, %d', chan);
for i = 3:length(s)
    eqn = strcat(eqn, ', :');
end
eqn = strcat(eqn, '));');

eqn = inline(eqn);
data = eqn(data);

% Get our data into readable variables
data_size = size(data);
np = data_size(1);
nd = length(data_size)-1;
max_points = zeros(1, nd);
max_points(:) = data_size(2:end);

% Initialize the current point and number of indirect points.
indir_points = prod(max_points);
point_size = num2cell(max_points);

s = [num_windows-1 point_size{:}];
out = zeros(s);

s(1) = num_windows;
points = zeros(s);
times = zeros(1, num_windows);

% Initialize the vectors to pull out the relevant information
% Each point in ms is 1000/sr ms.
% We want to start start ms in, so:
start_point = round(start*sr/1000);

% Now we need to know how often to sample. This is how many ms between
% pulses in ms. So it's calculated the same way as above.
inc_points = round(len*sr/1000);

% Finally we need to make a matrix containing vectors of points we want.
% Start with the window size in points:
window_size = round((sr/1000)*window);

% If the window size is bigger than the increment, that's an issue, so freak out.
if window_size > inc_points
    fprintf('Window is larger than the period! Not allowed.\n')
    return
end

% Now we calculate the margins on either side in points
margins = round((inc_points-window_size)*asym);

% Finally initialize the matrix, which is of size num_windows*2
% It will be the same for all dimension points, so this is just a 2 x num_windows matrix
% of spans.
windows = zeros(2, num_windows);
for i = 1:num_windows
    % First the start of the span
    windows(1, i) = start_point+inc_points*(i-1)+margins;

    % Then the end
    windows(2, i) = start_point+inc_points*(i-1)+margins+window_size;
end

% Make a plot so that you can see where you sampled from.
spans = zeros(np, 1);
for i = 1:length(windows)
    for j = windows(1, i):windows(2, i)
        if mod(i, 2) == 1
            spans(j) = 1;
        else
            spans(j) = -1;
        end
    end
end

% Generate the output vector.
for i = 1:indir_points
    [cp{1:length(max_points)}] = ind2sub(max_points, i);    
    % Get all the points.    
    for j = 1:num_windows
        points(j, cp{:}) = mean(data(windows(1, j):windows(2, j), cp{:}), 1);
    end 
    
    % Generate the output vector
    out(1:2:end, cp{:}) = points(1:2:end-1, cp{:})-points(2:2:end, cp{:});
    out(2:2:end, cp{:}) = -(points(2:2:end-1, cp{:})-points(3:2:end, cp{:}));
%     for j = 1:num_windows-1
%         out(j, cp{:}) = points(j+1, cp{:})-points(j, cp{:});
%     end
end

% Now generate the times vectors. These should correspond to the center
% point of "points". So basically the mean of the first dimension of
% windows.

times = mean(windows, 1);
	
	