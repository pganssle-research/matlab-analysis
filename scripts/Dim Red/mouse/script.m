function plot_mapping(input, map, discard_outliers)
    % Plots the mapping.

    if(nargin < 3)
        discard_outliers = 1;
    end
    
    vav = zeros(6, size(input, 2));
    col_map = [-2 -1 1 2 3 4];

    if(discard_outliers)
        for i = 1:6
            if(~isfield(map, 'conn_comp'))
                vav(i, :) = mean(input(cols(:) == col_map(i), :));
            else
                vav(i, :) = mean(input(cols(map.conn_comp) == col_map(i), :));
            end
        end
    else

        %Same thing, but discard outliers.
        noutliers = zeros(size(input, 1), 1);
        for i = 1:6
            if(~isfield(map, 'conn_comp'))
                lvec = cols(:) == col_map(i);
                ivec = input(lvec, :);
            else
                lvec = cols(map.conn_comp) == col_map(i);
                ivec = input(lvec, :);
            end

            mu = mean(ivec);
            dvec = arrayfun(@(x)(norm(mu-ivec(x, :))), 1:size(ivec, 1)); % Norm distance from mean
            sdev = std(dvec);

            not_outlier = abs(dvec) < 100*sdev;
            vav(i, :) = mean(ivec(not_outlier, :));

            noutliers(lvec) = not_outlier;

            if(sum(~not_outlier) > 0)
               fprintf('Set %d: Discarded %d outliers from %d points, %0.1f%% of data.\n', i, sum(~not_outlier), length(not_outlier), 100*sum(~not_outlier)/length(not_outlier));
            end   
        end
    end
    
    if(~exist('cm', 'var') && exist('cm.mat', 'file'))
        load('cm.mat', 'cm');
    end
    
    if(exist('cm', 'var'))
        figure(1);
        colormap(cm);
        figure(2)
        colormap(cm);
    end
    
    noutliers = logical(noutliers);
    if(size(vav, 2) > 2)
        figure(1);
        scatter3(input(noutliers, 1), input(noutliers, 2), input(logical(noutliers), 3), 50, cols(logical(noutliers)));
        figure(2);
        scatter3(vav(:, 1), vav(:, 2), vav(:, 3), 50, col_map);
    elseif(size(vav, 2) < 2)
        figure(1);
        scatter(input(logical(noutliers), 1), zeros(sum(noutliers), 1), 50, cols(logical(noutliers)));
        figure(2);
        scatter(vav(:), zeros(size(vav)), 50, col_map);
    else
        figure(1);
        scatter(input(logical(noutliers), 1), input(logical(noutliers), 2), 50, cols(logical(noutliers)));
        figure(2);
        scatter(vav(:, 1), vav(:, 2), 50, col_map);
    end
