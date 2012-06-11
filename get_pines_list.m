function [email, name, cat, url] = get_pines_list(tableID, html)

pageString = html;

% Pattern for finding MATLAB hooks
pattern = ['<a href="matlab:getTableFromWeb\(' num2str(tableID) '\)'];

% Find data from the table
[~, ~, ~, match] = regexp(pageString, [pattern '.*?<table[^>]*>.*?(<tr.*?>).*?</table[^>]*>' ], 'once');

emailn = {};
emaild = {};
name = {};
cat = {};
url = {};

m = 'notempty';
mmatch = match;
i = 1;
while(~isempty(m))
   [~, e, ~, m] = regexp(mmatch, '<tr.*?tr>', 'once');
   
   mmatch(1:e) = [];
      
   % Email address and URL
   nem = regexp(m, '<a href=\"javascript:reantispam\(''(?<first>.*?)''\s*,\s*''(?<last>.*?)''\);\">(?<name>.*?)</a>', 'names');
   urlm = regexp(m, '<a href=["''](?<url>http://.*?)["'']>.*?</a>', 'names');
   
   if(isempty(nem))
      if(~isempty(urlm))
         nem = regexp(m, '<td>\s*(?<name>.*?)\s*(\()?<a', 'names');
      else
         nem = regexp(m, '<td>\s*(?<name>.*?)\s*</td>', 'names'); 
      end
   end
   
   if(isempty(nem) || isempty(deblank([nem.name])))
      continue; 
   end
   
   u = [];
   cval = [];
   
   if(~isempty(urlm))
       u = urlm.url;
   end
   
   if(isfield(nem, 'first'))
        emailn{i} = deblank(nem.first);
        emaild{i} = deblank(nem.last);
   end
   
   url{i} = u;
   name{i} = nem.name;
   
   % Currently At
   catm = regexp(m, '<[Bb]>Currently At</[Bb]>:\s*(?<cat>.*?)\s*</div>', 'names');
   if(~isempty(catm))
        cval = catm.cat;
   end
   
   
   cat{i} = cval;
   
   i = i+1;
   
   if(i == 30)
      a = 4; 
   end
end

email = {emailn{:}; emaild{:}}';
