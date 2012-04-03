function plot_mc_struct(struct, point, fft, spans)
% Plots the struct array.
%
% Points will be interpreted as:
% [t, {dim_step}]
%
% Passing 0 to any entry will replace that with :
%   Leaving entries blank also replaces them with :, so 
% plot_mc_struct(struct, 0, 1) will plot all transients of all the FFTs.
% Passing -1 to an entry averages over that dimension (including transients)
% So plot_mc_struct(struct, [-1, 1]) plots the average of the first point
% in the first dimension (plus everything of the form out.adata(:, 1, :).
% 
% Pass -2 to plot the very first point in the remainder of the spectrum.
% So for input size (8000, 8, 25, 30, 30, 50), you can pass 
% [-1, 3, -2] to plot mean(in(:, :, 3, 1, 1, 1), 2).
%
% Pass -3 to a point to plot all the rest of the spectrum, so for
% [-1, 3, 2, -3], you'd plot mean(in(:, :, 3, 2, :), 2)
%
% Default is [-1, -2], pass anything below -2 or anything that's not an int
% to for this behavior.

if(~exist('point', 'var') || ~isnumeric(point) || ~isempty(find(point < -3, 1, 'first')))
   point = [-1, -2]; 
end

e = logical(arrayfun(@(x)(x == -2 || x == -3), point));

l = find(e, 1, 'first');
if(~isempty(l) && l < length(point))
    point((l+1):end) = [];
end

if(~exist('fft', 'var') || ~fft)        
    data = struct.mdata;
    x = struct.t;
else
    if(~isfield(struct, 'fft'))
        struct = add_fft(struct);
    end
    
    data = struct.fft;
    x = struct.f;
end

for i = find(point == -1)
    data = mean(data, i+1);
end 

cmd = 'squeeze(data(:';
j = 0;
l = length(size(data));
for i = 1:length(point)
    if(point(i) > 0)
        cmd = sprintf('%s,%d', cmd, point(i));
        j = j+1;
    elseif(point(i) > -2)
        cmd = [cmd ',:'];  %#ok
        j = j+1;
    elseif(point(i) == -2)
        cmd = [cmd, repmat(sprintf(',%d', 1), 1, l-i)];  %#ok
        j = j+length(size(data))-i;
        break;
    else
        break;
    end
end

if(j < l)
    cmd = [cmd, repmat(',:', 1, l-j-1)];
end

cmd = [cmd, '))'];
data = eval(cmd);

if(fft == 2)
    data = abs(data);
end

if(~fft && exist('spans', 'var'))
   if(size(spans, 1) ~= size(data, 1))
       if(size(spans, 2) == size(data, 1))
           spans = spans';
       else
           clear spans;
       end
   end
   
   if(exist('spans', 'var'))
    data = [data(:, :), spans]; 
   end
end

plot(x, data(:, :));









