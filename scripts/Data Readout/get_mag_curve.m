function [out, spans, points, times] = get_mag_curve(curve, len, num_windows, sr, window, data, asym, chan)
% Function for getting the magnetization curve as a function of time from an
% experiment consisting of a train of 180s.
%
% start = Time point to start from (ms)
% len = length of the plateaus (ms)
% num_windows = number of 'points' to acquire
% window = how long to average between pulses (ms)
% sr is the sampling rate -> If this is set to <0 and data is a struct,
% takes data.prog.sr
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
% [out, spans, points, times] = get_mag(start, len, num_windows, sr, window, data, asym, chan);
% out = get_mag(curve); -> If you already have the curve

if(nargin > 1) 
    if(isstruct(data) && sr <= 0)
        sr = data.prog.sr;
    end
    
    if(nargin < 6)
        error('Insufficient number of points.');
    end
    
    start = curve;
    
    if(nargin < 7)
        [curve, spans, points, times] = get_mag(start, len, num_windows, sr, window, data);
    elseif(nargin < 8)
        [curve, spans, points, times] = get_mag(start, len, num_windows, sr, window, data, asym);
    else
        [curve, spans, points, times] = get_mag(start, len, num_windows, sr, window, data, asym, chan);
    end
else
    spans = 0;
    points = 0;
    times = 0;
end

out = curve(1:2:end-1, :)-curve(2:2:end, :);





