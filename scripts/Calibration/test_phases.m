function test_phases(in, freq, start_phase, end_phase, num)
% Generates a new bunch of plots with various phases  for phasing
% spectra
%
% in = Complex data
% freq = frequency vector (defaults to 1:np)
% start_phase = first phase you want to look at (defaults to 0)
% end_phase = last phase you want to see (defaults to 2pi)
% num = number of plots you want (defaults to 16)
%
% Usage:
% test_phases(in, freq, start_phase, end_phase, num);

% Set up the defaults
if nargin < 2
    freq = 1:np;
end

if nargin < 3 || start_phase <= 0
    start_phase = 0;
end

if nargin < 4 || end_phase <= 0
    end_phase = 2*pi;
end

if nargin < 5 || num <= 1;
    num = 16;
end

if end_phase > 2*pi
    end_phase = rem(end_phase, 2*pi);
end

if start_phase > 2*pi
    start_phase = rem(start_phase, 2*pi);
end

step_size = (end_phase - start_phase)/(num-1);
phases = start_phase:step_size:end_phase;
phase_vec = transpose(exp(1i*phases)); % Generate a column vector with the phases we want to try

% Make sure the data vector is a 1D row vector
s = size(in);

if(numel(s) > 2)
    fprintf('Must be 1D vector.\n');
end

% Assuming input data is direct x transients
data = transpose(mean(in, 2));

% Do all the plotting
new_vectors = real(phase_vec*data);
figure
handles = plot(freq, new_vectors);
xlabel('Frequency (\omega)');


for j = 1:numel(handles)
    set(handles(j), 'DisplayName', ['Phase = ', num2str(phases(j))]);
end
