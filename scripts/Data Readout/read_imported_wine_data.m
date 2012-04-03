% Load up all the wine data in a folder.
function [all_data, data_1d, data_2d] = read_imported_wine_data(path)

if(~exist('path', 'var') || ~exist(path, 'file'))
    path = uigetdir2('C:\Users\Omega\Dropbox\');
end

if(path(end) ~= '\')
    path = strcat(path, '\');
end

listing = dir(path);
listing(1:2) = [];

j = 1;

proto_struct = struct('Title', '', 'is2d', logical(0),  'Info', struct(), 'Acqus', struct(), ...
    'Data', [], 'IData', [], 'XAxis', [], 'YAxis', [], ...
    'Proc2s', struct(), 'Acqu2s', struct());

names = fieldnames(proto_struct);

for i = 1:length(listing)
    l = listing(i);
    if(l.isdir)
        continue;
    end

    load(strcat(path, l.name), '-mat', 'd');

    if(~exist('d'))
        continue;
    end

    all_data(j) = proto_struct;
        
    for k = 1:length(names)
        n = char(names(k));
       if(eval(['isfield(d, ''' n ''')']))
          eval(['all_data(j).' n ' = d.' n ';']);
       end
    end
    
    if(isfield(d, 'YAxis'))
       all_data(j).is2d = logical(1); 
    end
    
    j = j+1;

    clear d;   
end

% Do a little data processing.
data_1d = all_data(~[all_data.is2d]);
data_2d = all_data([all_data.is2d]);
    