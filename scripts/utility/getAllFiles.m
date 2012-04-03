function [fileList, fullPath] = getAllFiles(dirName)
  
  if(dirName(end) ~= '\')
      dirName = strcat(dirName, '\');
  end
  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  if(~isempty(fileList))
    fullPath = strcat(dirName, fileList);
  else
      fullPath = fileList;
  end
  
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    [fl, fp] = getAllFiles(nextDir);    %# Recursively call getAllFiles]
    fileList = [fileList; fl];
    fullPath = [fullPath; fp];
  end
  
  

end