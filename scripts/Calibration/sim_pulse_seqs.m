function out = sim_pulse_seqs(init, fin, num_steps, pi_val, num_pulses) 
% Generates a simulation of the output of an experiment where a simple
% 180-wait-180 experiment is performed while varying the 180 time, for
% comparison to real values.
%
% init_pulse = time in us, default = 0.0
% fin_pulse = time in us, default = 250.0
% num_steps = number of steps, default = 30
% pi_val = true pi time, default = 120.0
% num_pulses = how many pulses to show, default = 50

if(~exist('init', 'var'))
   init = 0.0;
end

if(~exist('fin', 'var'))
   fin = 250.0; 
end

if(~exist('num_steps', 'var') || num_steps < 2)
    num_steps = 30;
end

if(~exist('pi_val', 'var') || pi_val <= 0.0)
    pi_val = 120.0;
end

if(~exist('num_pulses', 'var') || num_pulses < 1)
    num_pulses = 50;
end

pulse = pi*linspace(init, fin, num_steps)/pi_val; % As a percentage of a pi pulse
steps = reshape(repmat(1:num_pulses, 60, 1), 60*num_pulses, 1);

out = cell2mat(arrayfun(@(x)cos(pulse.*x), steps, 'UniformOutput', false));

end