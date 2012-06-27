function plot_mc_struct(struct, point, fft, spans, first, span)
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
%
% If 'first' is a scalar logical positive, it interprets all positive 
% values as, "plot the first (x) entries in this dimensions. If it's a 
% vector, it's element-wise evaluation.
%
% Usage:
% plot_mc_struct(struct, point, fft, spans, first);

if(~exist('point', 'var') || ~isnumeric(point) || ~isempty(find(point < -3, 1, 'first')))
   point = [-1, -2]; 
end

if(~exist('first', 'var'))
    first = false;
else
    first = logical(first);
end

if(isscalar(first))
   first = repmat(first, size(point));
elseif(length(first) < length(point))
   first(end:length(point)) = false; 
end

e = logical(arrayfun(@(x)(x == -2 || x == -3), point));

l = find(e, 1, 'first');
if(~isempty(l) && l < length(point))
    point((l+1):end) = [];
end

if(~exist('fft', 'var') || ~fft)        
    data = struct.mdata;
    x = struct.t;
	 xlab = 'Time (s)';
else
    if(~isfield(struct, 'fft'))
        struct = add_fft(struct);
    end
    
    x = struct.f; 
	 xlab = 'Frequency (Hz)';
    
    data = struct.fft;
    
    if(fft == 1)
        if(isfield(struct, 'disp') && isfield(struct.disp, 'phase') && length(struct.disp.phase) == 3)
            phase = fliplr(pi*struct.disp.phase/180);
            n = size(data(:, :), 2);
            phasecor = (exp(-1i*polyval(phase, x)))';

            for i = 1:n
                data(:, i) = data(:, i).*phasecor;
            end
        
        end
        
        data = real(data);
    end
end

for i = find(point == -1)
    data = mean(data, i+1);
end 
cmd = 'squeeze(data(:';

j = 0;
l = length(size(data));
for i = 1:length(point)
    if(point(i) > 0)
        if(~first(i))
            cmd = sprintf('%s,%d', cmd, point(i));
        else
            cmd = sprintf('%s,1:%d', cmd, point(i));
        end
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

ylab = 'Signal (V)';
if(isfield(struct, 'disp') && isfield(struct.disp, 'mag_cal'))
	data = data*struct.disp.mag_cal;
	ylab = 'Magnetic field (pT)';
end

plot(x, data(:, :));

xlabel(xlab);
ylabel(ylab);

if(isfield(struct, 'FileName'))
   title(struct.FileName); 
end

if(exist('span', 'var'))
    if(span(1) > span(2))
        xmax = span(1);
        xmin = span(2);        
    elseif(span(2) > span(1))
        xmax = span(2);
        xmin = span(1);        
    else
        return;
    end
    
    [~, mini] = min(abs(x-xmin));
    [~, maxi] = min(abs(x-xmax));
    
    if(mini == maxi)
        return;
    end
    
    ymax = max(max(data(mini:maxi, :)));
    ymin = min(min(data(mini:maxi, :)));
    
    padding = (ymax-ymin)*0.05;
    
    axis([xmin, xmax, ymin-padding, ymax+padding]);
end