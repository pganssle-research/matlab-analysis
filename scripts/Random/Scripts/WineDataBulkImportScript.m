% Bulk import of Wine Data.

path = uigetdir2('C:\Users\Omega\Dropbox\');
listing = dir(path);
listing(1:2) = [];

if(path(end) ~= '\') 
   path = strcat(path, '\');
end

datapath = strcat(path, 'MatlabData\');

if(~exist(datapath, 'file')) 
   mkdir(datapath); 
end

j = 1;
for i = 1:length(listing)
    if(~listing(i).isdir)
        continue;
    end
   
    name = strcat(path, listing(i).name);
   
    % Search for 1r or 2rr.
    names = getAllFiles(name);
    if isempty(find(strcmp(names, '1r'), 1)) && isempty(find(strcmp(names, '2rr'), 1))
        continue;
    end  
    
    d = rbnmr(name);

    j = j + 1;
    
    matname = sprintf('%s%s.mat', datapath, listing(i).name);
    save(matname, 'd');    
end

all_data = read_imported_wine_data(datapath);
