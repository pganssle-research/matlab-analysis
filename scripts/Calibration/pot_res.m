function out = pot_res(tot, pos, para)
% Gives potentiometer resistance at position pos (0-100) with optional
% resistor para in parallel.
%
% out = pot_res(tot, pos, para);

pos = pos/100;
out = pos.*tot;

if(exist('para', 'var'))
    out = (out*para)./(out+para);
end

