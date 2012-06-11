newText = getTableFromURLs('http://waugh.cchem.berkeley.edu/global_nuts.html');

[email, name, cat, url] = get_pines_list(2, newText);

% Generate the file
nl = sprintf('\n');

file = ['<?xml version="1.0"?>' nl '<nuts>' nl];


for i = 1:size(name, 2)
    nstring = [];
    
    name{i} = strrep(name{i}, '  ', ' ');
    name{i} = strrep(name{i}, '&nbsp;', '');
    
    nstring = ['<pinenut' nl 'name="' name{i} '"' nl];
    
    if(~isempty(email{i, 1}) && ~isempty(email{i, 2}))
        nstring = [nstring 'emailn="' email{i, 1} '" emaild="' email{i, 2} '"' nl];
    end
    
    if(~isempty(cat{i}))
       nstring = [nstring 'cat="' cat{i} '"' nl]; 
    end
    
    if(~isempty(url{i}))
       nstring = [nstring 'url="' url{i} '"' nl]; 
    end
    
    nstring = [nstring '/>' nl];
    
    file = [file nstring nl];    
end

file = [file '</nuts>'];

fid = fopen('global_nuts.xml', 'w');
fprintf(fid, '%s', file);
fclose(fid);
