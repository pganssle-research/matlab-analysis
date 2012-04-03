% Generate random data.
total_dsets = 20;
dsets_1 = 6;

t = 0:0.001:10-0.001;
dsets = zeros(total_dsets, length(t));

% These are lorentzian weighted
for j=1:dsets_1
    for k = 1:2
        dsets(j, :) = dsets(j, :) + exp(-(t/(5+rand(1)))).*exp(1i*2*450*pi*rand(1)*t);
    end
    dsets(j, :) = dsets(j,:)/sum(dsets(j, :));
end

% These are Gaussian weighted
for j=(dsets_1+1):total_dsets
    for k=1:10
        dsets(j, :) = dsets(j, :) + exp(-(t/(rand(1)))).*exp(1i*2*(450*rand(1))*pi*t);
    end
    dsets(j, :) = dsets(j,:)/sum(dsets(j, :));
end

% Fourier transform
sr = 1/(t(2)-t(1));
f = 0:(sr/(length(t)-1)):(sr/2);
sdsets = fft(dsets, [], 2);

% Test a bunch of these maps
techs = {'PCA', 'LDA', 'MDS', 'FactorAnalysis', 'GPLVM', 'Sammon', 'IsoMap'};
mapped_sdset = compute_mapping(sdsets, 'LDA');
mapped_dset = compute_mapping(dsets, 'LDA');

scatter(mapped_sdset(:, 1), mapped_sdset(:, 2), 'MarkerEdgeColor', 'r', 'Title', 'BlahBlah')
