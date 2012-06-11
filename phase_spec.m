figure
span = [130, 150];
hold all

num = 16;
start = 0;
fin = 360;

out.disp.phase = [0, 0, 0];
phases = linspace(start, fin, num);
strings = num2cell(phases);

for i = 1:num
   out.disp.phase(1) = phases(i);
   plot_mc_struct(out, [-2], 1, -1, -1, span);
   strings{i} = sprintf('%f', strings{i});
end

legend(strings{:});

clear span;
clear num
clear start
clear fin
clear strings;