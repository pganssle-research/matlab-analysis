function [out, stdev] = make_avg_anything(data, end_point, start_point)
% Takes the average ratio across the curve and scales, then averages them.
%
% Usage: [out, stdev] = make_avg_anything(data, end_meas, start_meas, end_point, start_point);
%
% all_flag: If 1, include all transients in the average, even if they aren't included in making the average ratio. Default to 0.
% end_point: The last point to be included in the averaging, defaults to the last point (pass 0 omit for default)
% start_point: The first point to be included, defaults to first point (pass 0 or omit for default)
%
% Example, with two curves a1 = [0.5 0.3 0.2 0.1 0.05], a2 = [0.24 0.155 0.105 0.053 0.026]
% It will generate the scale vector [2.0833 1.9355 1.9048 1.8868 1.9231]
% The mean of this is 1.9467, so it will scale everything in the second vector up by 1.9467:
%
% a2 = [0.4672 0.3017 0.2044 0.1032 0.0506]
%
% Then it will average a1 and a2 and give standard deviations:
%
% out = [0.4836 0.3009 0.2022 0.1016 0.0503]
% stdev = [0.0232 0.0012 0.0031 0.0023 0.0004]

s = size(data);

% Set to defaults if necessary
if nargin < 2 || end_point < 1 || end_point > s(1)
	end_point = s(1);
end

if nargin < 3 || start_point < 1 || start_point > end_point
	start_point = 1;
end

% Scale all the stuff appropriately
s = num2cell(s);

if(size(s) > 1)
    out = zeros(s{2:end});
    stdev = zeros(s{2:end});
else
    out = zeros(s{2:end}, 1);
    stdev = zeros(s{2:end}, 1);
end

out(:) = arrayfun(@(x)mean(data(start_point:end_point, x)), 1:length(data(1, :)));
stdev(:) = arrayfun(@(x)mean(data(start_point:end_point, x)), 1:length(data(1, :)));
% 
% new_data(:, 1) = data(:, 1);
% for i=2:s(2)
% 	avg_vec = data(start_point:end_point, 1)./data(start_point:end_point, i);
% 	avg = mean(avg_vec);
% 	
% 	new_data(:, i) = data(:, i)*avg;
% end
% 
% % Generate the outputs now
% out = zeros(s(1), 1);
% stdev = zeros(s(1), 1);
% 
% for i=1:s(1)
%     out(i) = mean(new_data(i, :));
%     stdev(i) = std(new_data(i, :));
% end

