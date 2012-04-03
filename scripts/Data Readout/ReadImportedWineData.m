% Load up all the wine data in a folder.
function all_data = read_imported_wine_data(path)
if(~exist('path', 'var') || ~exist(path, 'file'))
    path = uigetdir2('C:\Users\Omega\Dropbox\');
end

listing = dir(path);
listing(1:2) = [];