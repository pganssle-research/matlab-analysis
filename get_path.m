function path = get_path(histpath, filetypes, last)
% Simple way to generate UI popups to get paths, with history.
%
% Default histpath is get_path_hist, generic.
% Default filetypes is {'*.*'} - should be a cell delimited by ;
%
% Usage:
% path = get_path([histpath, filetypes]);

if(~exist('histpath', 'var') || isempty(histpath) || ~ischar(histpath))
    histpath = 'get_path_hist.mat';
end

if(~exist('last', 'var'))
    last = 0;
end

history = {};
phistory = {};

if ~exist(histpath, 'file')
    default_dir = 'C:\Omega\Data\'; % If the readout_history file is missing,
    if(~exist(default_dir, 'file'))
        default_dir = pwd;
    end
else
    load(histpath);
    
    default_dir = pwd;
    
    
    if (~isempty(history))
        if(last < 0 && ~isempty(phistory))
           bitset = 0;
           
           last = abs(last);
           
           k = 1;
           
           for i = length(phistory):-1:1
               if(exist(phistory{i}, 'file'))
                   bitset = 1;
                   path = phistory{i};
                    
                   k = k+1;
                   
                   if(k >= last)
                      break; 
                   end
               end
           end
            
           if(~bitset)
               last = 0;
           else
               lfs = find(path == filesep, 1, 'last');
               if(~isempty(lfs))
                  filefolder = path;
                  filefolder((lfs+1):end) = [];
                  
                  update_history(histpath, path, filefolder, history, phistory);
               end
               
               return;
           end
        end
        
        % Now search through the history file to find the most
        % recently used one that exists. This is useful if the same
        % history file is synced across multiple systems. We'll set
        % a variable keepatmost in the readout_history.mat file, so
        % that we can adjust how long the history we want to keep
        % is. Default is keep all., keepatmost == -1 also means
        % keep all.
        k = 1;
        for j = length(history):-1:1
            if exist(history{j}, 'file')
                default_dir = history{j};
                k = k + 1;

                if(k >= last)
                    break; % Stop looking once you've found the 'last'-th one.
                end
            end
        end
        dupes = ismember(history, default_dir); % List of the positions of duplicate entries
        dupes = dupes(1:end-1); % The most recent one is OK to stay, the others shouldn't even be there.

        history(dupes) = []; %#ok<*NODEF> Delete the relevant entries
    end
end

if(~exist('filetypes', 'var') || ~iscell(filetypes))
    filetypes = {'*.*'};
end

[filepath,filefolder]=uigetfile(filetypes,'Select a file', default_dir);
path=fullfile(filefolder,filepath);

if (isempty(path) || (~exist(path, 'file')))
    warning('Invalid file name.'); %#ok
    return;
end

update_history(histpath, path, filefolder, history, phistory);


function update_history(histpath, path, filefolder, history, phistory)
% Update the history
if(~isempty(filefolder) && ischar(filefolder))
    history(ismember(history, filefolder)) = [];
    history = [history, {filefolder}];
end

if(~isempty(path) && ischar(path))
    phistory(ismember(phistory, path)) = [];
    phistory = [phistory, {path}];
end

if(exist('keepatmost', 'var') && keepatmost >= 1)
    if(length(history) > keepatmost)
        history = history((end-keepatmost):end);
    end
    
    if(length(phistory) > keepatmost)
        phistory = phistory((end-keepatmost):end);
    end
end

varlist = {};

if(~isempty(history))
    varlist = [varlist, {'history'}];
end

if(~isempty(phistory))
    varlist = [varlist, {'phistory'}];
end

if(exist(histpath, 'file'))
    save(histpath, varlist{:}, '-append');
else
    save(histpath, varlist{:});
end