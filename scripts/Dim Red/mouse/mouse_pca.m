% Run all_mouse_data first.

n = 189;
% Now we need to generate an array of the subsection from ppm = -0.5 to ppm = -5.
[v1, ind1] = min(abs(mice(1).p+0.5));
[v2, ind2] = min(abs(mice(1).p+4.5));
np = (ind2-ind1)+1;

% Rows are observations, columns are dimensions.
m_pca_inr = zeros(n, np);
m_pca_ini = zeros(n, np);
m_pca_inr(1, :) = real((mice(1).fft(ind1:ind2)))';
m_pca_ini(1, :) = imag((mice(1).fft(ind1:ind2)))';

for i = 2:n
%     [~, ind1] = min(abs(mice(i).p+0.5));
%     [~, ind2] = min(abs(mice(i).p+5.5));
%     
     if((ind2-ind1)+1 ~= np)
         dp = (ind2-ind1+1)-np;
         disp(sprintf('Disparity, size %d\n', dp));
         ind1 = ind1 + dp;
     end
    
    m_pca_inr(i, :) = real(mice(i).fft(ind1:ind2))';
    m_pca_ini(i, :) = imag(mice(i).fft(ind1:ind2))';
end

% Center data, then scale to unit variance.
m_pca_in = zeros(n, np);
m_pca_in(:, :) = magnitude(m_pca_inr+m_pca_ini);
%m_pca_in(:, (np+1):end) = m_pca_ini;

m_pca_in = bsxfun(@minus, m_pca_in, mean(m_pca_in));    % Zero mean
%m_pca_in = bsxfun(@rdivide, m_pca_in, std(m_pca_in));   % Unit variance

[pc_o pmap] = compute_mapping(m_pca_in, 'PCA', 0.9);
[nc_o nmap] = compute_mapping(m_pca_in, 'HessianLLE', 2, 20, 'JDQR');
[im_o imap] = compute_mapping(m_pca_in, 'Isomap', 3, 30);

cols = zeros(1, n);
cols(find(strncmp('spec_A', {mice.name}, 6))) = 2;
cols(find(strncmp('spec_B', {mice.name}, 6))) = 3;
cols(find(strncmp('spec_C', {mice.name}, 6))) = -1;
cols(find(strncmp('spec_D', {mice.name}, 6))) = 1;
cols(find(strncmp('spec_E', {mice.name}, 6))) = 4;
cols(find(strncmp('WT', {mice.name}, 2))) = -2;

% Save this shit.
save('mice.mat', 'mice', 'm_pca_in', 'pc_o', 'pmap', 'cols');