%% Get data from an HTML table straight into MATLAB
% Copyright 2008 - 2010 The MathWorks, Inc.

%% Navigate to a webpage with a table 
% You have to use the MATLAB browser because this takes advantage of the
% matlabcolon (matlab:) protocol.
artists = {};
titles = {};
length = {};

albums = 2:41;
albums(albums==24) = [];

for i = albums
    url = sprintf('http://en.wikipedia.org/wiki/Now_That\''s_What_I_Call_Music!_%d_(U.S._series)', i)
    %% Call the function to grab the data
    newText = getTableFromURLs(url);
    myTableData = [];
    
    if(i == 26)
        k = 2;
    else
        k = 3;
    end
    
    while(size(myTableData, 2) ~= 4)
        myTableData = getTableFromURLs(newText, k);
        k = k+1;
    end
    %% Updated Browser
    % This shows what the browser looks like.  This is just a screenshot but it
    % would be whatever page you are on with MATLAB icons next to tables.
    %% Select the table 
    % You can select the table interactively by clicking on the MATLAB icon or
    % passing a number to the getTableFromWeb function.  I have selected to get
    % the Valuation Measures
    
    now3 = cell2struct(myTableData(2:end, 2:end)', {'Artist', 'Title', 'Length'}, 1);
    
    artists = {artists{:}, now3(:).Artist};
    titles = {titles{:}, now3(:).Title};
    length = {length{:}, now3(:).Length};
    now = struct('Artist', artists, 'Title', titles, 'Length', length);
end

