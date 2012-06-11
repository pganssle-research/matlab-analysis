function Parameters = ReadParameterFile2(FileName)
% Generates a cell array of all the lines in the procpar file.

fp = fopen(FileName, 'r');
if (fp == -1)
  disp(['Couldn''t open parameter file "' FileName '". Aborting ...']);
  Parameters = [];
  return
end

%
%read first line and edit
%
tmp = fgetl(fp);
tmp([strfind(tmp, '#') strfind(tmp, '$') strfind(tmp, ';')]) = '';
Parameters = {tmp};

while 1
  tmp = fgetl(fp);
  if ischar(tmp)
    tmp = deblank(tmp);
    if (length(tmp)<80)		%avoid endless lines and thus and memory time loss
      tmp([strfind(tmp, '#') strfind(tmp, '$')]) = '';
      Parameters = [Parameters, {tmp}];
    end
  else
    break
  end
end

fclose(fp);
