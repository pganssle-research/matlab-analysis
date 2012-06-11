function out = adjust_baseline(in, signal_loc)
% Give this an mc_struct and it will alter the data to subtract off a
% running-average of the signal. If signal_loc is provided, the window will
% be 10x wider than a single period. Default signal_loc is 100Hz.
%
% signal_loc must be in Hz.
%
% Usage:
% out = adjust_baseline(in[, signal_loc]);

out = in;

if(~exist('signal_loc', 'var'))
    signal_loc = 100;
end

windowSize = round((in.prog.sr/signal_loc)*10);   % Window size in points;
perSide = round(windowSize/2);

if(~isfield(out, 'odata'))
    out.odata = out.mdata;
end


done = perSide;
if(perSide > in.prog.np)
    for i = 1:in.prog.np
        out.mdata(i, :) = in.mdata(i, :) - mean(in.mdata(:, :), 1);
    end
elseif(windowSize > in.prog.np)
    done = in.prog.np-perSide;
end

for i = 1:done
   out.mdata(i, :) = in.mdata(i, :) - mean(in.mdata(1:(i+perSide), :), 1);
end

for i = ((done+1):(in.prog.np-perSide))
    out.mdata(i, :) = in.mdata(i, :) - mean(in.mdata((i-perSide):(i+perSide), :));
end

for i = ((in.prog.np-perSide-1):in.prog.np)
   out.mdata(i, :) = in.mdata(i, :) - mean(in.mdata((i-perSide):end, :)); 
end


