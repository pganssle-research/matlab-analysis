function out = apodize1d(data, apo, sr)
% Apodizes a set of 1D spectra by applying an exponential apodization to
% each individual spectrum.
%
% Usage:
% out = apodize1d(data, apo, sr);

if(apo <= 0)
    out = data;
    return;
end

s = num2cell(size(data)); % Preserve size
t = (1:s{1})'/sr; % Time vector
out = cell2mat(arrayfun(@(x)data(:, x).*exp(-t/apo), 1:length(data(1, :)), 'UniformOutput', false)');
out = reshape(out, s{:});