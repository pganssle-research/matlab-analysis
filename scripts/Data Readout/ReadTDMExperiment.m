function struct = ReadTDMExperiment(filename)

%Recreate needed property constants defined in nilibddc_m.h
DDC_FILE_NAME					=	'name';
DDC_FILE_DESCRIPTION			=	'description';
DDC_FILE_TITLE					=	'title';
DDC_FILE_AUTHOR					=	'author';
DDC_FILE_DATETIME				=	'datetime';
DDC_CHANNELGROUP_NAME			=	'name';
DDC_CHANNELGROUP_DESCRIPTION	=	'description';
DDC_CHANNEL_NAME				=	'name';

% Use the property constants defined in the Magnetometer Controller program
MCTD_MAINDATA = 'DataGroup';            % The main data group name.
MCTD_AVGDATA = 'AvgGroup';				% The group of data averages.
MCTD_PGNAME = 'ProgramGroup';             % Name of the program group

MCTD_NP = 'NPoints';                    % Number of points (data)
MCTD_SR = 'SamplingRate';				% Sampling rate 
MCTD_NT = 'NTransients';				% Number of transients total
MCTD_NDIM = 'NDims';					% Number of dimensions
MCTD_NCYCS = 'NCycles';                 % Number of cycles
MCTD_NCHANS = 'NChans';                 % Number of channels
MCTD_TIMESTART = 'TimeStarted';         % Time started
MCTD_TIMEDONE = 'TimeCompleted';		% Time completed

MCTD_CSTEP = 'CurrentStep';             % Current step as a linear index
MCTD_CSTEPSTR = 'CurrentStepString';	% Current step as a string

MCTD_PNP = 'nPoints';                   % p.np
MCTD_PSR = 'SampleRate';                % p.sr
MCTD_PNT = 'nTransients';               % p.nt

MCTD_PTMODE = 'TransAcqMode';			% p.tfirst
MCTD_PNINSTRS= 'nInstrs';				% p.n_inst
MCTD_PNUINSTRS= 'nUniqueInstrs';		% p.nUniqueInstrs

MCTD_PVAR= 'Varied';					% p.varied
MCTD_PNVAR= 'nVaried';                  % p.nVaried
MCTD_PNDIMS= 'nDims';					% p.nDims
MCTD_PNCYCS= 'nCycles';                 % p.nCycles
MCTD_PSKIP= 'Skip';                     % p.skip
MCTD_PMAX_N_STEPS= 'MaxNSteps';         % p.max_n_steps
MCTD_PREAL_N_STEPS= 'RealNSteps';		% p.real_n_steps

MCTD_NFUNCS= 'nFuncs';                  % p.nFuncs
MCTD_TFUNCS= 'tFuncs';                  % p.tFuncs

% Program Instruction Details
MCTD_PROGFLAG= 'ProgFlag';              % p.instrs[i].flags;
MCTD_PROGTIME= 'ProgTime';              % p.instrs[i].instr_delay
MCTD_PROGINSTR= 'ProgInstr';			% p.instrs[i].instr;
MCTD_PROGUS= 'ProgUS';                  % p.instrs[i].units and p.instrs[i].scan
MCTD_PROGID= 'ProgID';                  % p.instrs[i].instr_data;

% Varying program details
MCTD_PMAXSTEPS= 'MaxSteps';             % p.maxsteps
MCTD_PVINS= 'VarIns';					% p.v_ins
MCTD_PVINSDIM= 'VarInsDims';			% p.v_ins_dim
MCTD_PVINSMODE= 'VarInsMode';			% p.v_ins_mode
MCTD_PVINSLOCS= 'VarInsLocs';			% p.v_ins_locs
MCTD_PSKIPLOCS = 'SkipLocs';             % p.skip_locs
MCTD_PSKIPEXPR= 'SkipExpr';             % p.skip_expr
MCTD_PDELEXPRS= 'DelayExprs';			% p.delay_exprs
MCTD_PDATEXPRS= 'DataExprs';			% p.data_exprs

MCTD_PFUNCLOCS= 'FuncLocs';             % p.func_locs

% Function details
MCTD_PF_NAME= 'FuncNames';              % pfunc.name
MCTD_PF_RF= 'FuncResFlags';             % pfunc.r_flags
MCTD_PF_NINSTR= 'FuncNInstrs';          % pfunc.n_instr
MCTD_PF_ARGMODE= 'FuncArgMode';         % pfunc.argmode

% Function instruction details
MCTD_PFINSTR= 'FuncInstr';              % pfunc.instrs[i].instr
MCTD_PFFLAG= 'FuncFlag';				% pfunc.instrs[i].flags
MCTD_PFID= 'FuncID';					% pfunc.instrs[i].instr_data;
MCTD_PFUS= 'FuncUS';					% pfunc.instrs[i].time_units and scan

%Check if the paths to 'nilibddc.dll' and 'nilibddc_m.h' have been
%selected. If not, prompt the user to browse to each of the files.
hist = {};
NI_TDM_DLL_Path = 'C:\Users\Omega\Documents\MATLAB\dev\bin\64-bit\nilibddc.dll';
NI_TDM_H_Path = 'C:\Users\Omega\Documents\MATLAB\dev\include\64-bit\nilibddc_m.h';

if(exist('ReadTDMExperimentHistory.mat', 'file'))
    load('ReadTDMExperimentHistory.mat');
end

if exist(NI_TDM_DLL_Path,'file')==0
    [dllfile,dllfolder]=uigetfile('*dll','Select nilibddc.dll');
    NI_TDM_DLL_Path=fullfile(dllfolder,dllfile);
end

if exist(NI_TDM_H_Path,'file')==0
    [hfile,hfolder]=uigetfile('*h','Select nilibddc_m.h');
    NI_TDM_H_Path=fullfile(hfolder,hfile);
end

[~, libname, ~] = fileparts(NI_TDM_DLL_Path);

%Prompt the user to browse to the path of the TDM or TDMS file to read
if(nargin < 1 || ~exist(filename, 'file'))  
    if ~exist('ReadTDMExperimentHistory.mat', 'file')
        default_dir = 'C:\Omega\Data\'; % If the readout_history file is missing,
        if(~exist(default_dir, 'file'))
            default_dir = pwd;
        end
    else
        if (length(hist) < 1) 
            default_dir = pwd;
        else
            % Now search through the history file to find the most
            % recently used one that exists. This is useful if the same
            % history file is synced across multiple systems. We'll set
            % a variable keepatmost in the readout_history.mat file, so
            % that we can adjust how long the history we want to keep
            % is. Default is keep all., keepatmost == -1 also means
            % keep all.
            default_dir = pwd;
            for j = length(hist):-1:1 
                if exist(hist{j}, 'file')
                    default_dir = hist{j};
                    break; % Stop looking once you've found it
                end
            end
            dupes = ismember(hist, default_dir); % List of the positions of duplicate entries
            dupes = dupes(1:end-1); % The most recent one is OK to stay, the others shouldn't even be there.

            hist(dupes) = []; %#ok<*NODEF> Delete the relevant entries 
        end  
    end
    
    [filepath,filefolder]=uigetfile({'*.tdm';'*.tdms'},'Select a TDM or TDMS file', default_dir);
    Data_Path=fullfile(filefolder,filepath);
    struct.path = Data_Path;
    
    if (isempty(Data_Path) || (~exist(Data_Path, 'file')))
        warning('Invalid file name.'); %#ok
        return;
    end

    % Update the history
    hist(ismember(hist, filefolder)) = [];
    hist = [hist, {filefolder}];
    
    if(exist('keepatmost', 'var') && keepatmost >= 1 && length(hist) > keepatmost) 
       hist = hist((end-keepatmost):end); %#ok 
    end
   
    if(~isempty('hist'))
        if(exist('ReadTDMExperimentHistory.mat', 'file'))
            save('ReadTDMExperimentHistory.mat', 'hist', 'NI_TDM_DLL_Path', 'NI_TDM_H_Path', '-append');
        else
            save('ReadTDMExperimentHistory.mat', 'hist', 'NI_TDM_DLL_Path', 'NI_TDM_H_Path');
        end  
    end
else
    Data_Path=filename;
    struct.path = Data_Path;
end
%Load nilibddc.dll (Always call 'unloadlibrary(libname)' after finished using the library)
loadlibrary(NI_TDM_DLL_Path,NI_TDM_H_Path);

%Open the file (Always call 'DDC_CloseFile' when you are finished using a file)
fileIn = 0;
[err,dummyVar,dummyVar,file]=calllib(libname,'DDC_OpenFileEx',Data_Path,'',1,fileIn);

%Read and display file name property
filenamelenIn = 0;

%Get the length of the 'DDC_FILE_NAME' string property
[err,dummyVar,filenamelen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_NAME,filenamelenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfilename=libpointer('stringPtr',blanks(filenamelen));
    [err,dummyVar,filename]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_NAME,pfilename,filenamelen+1);
    setdatatype(filename,'int8Ptr',1,filenamelen);
    struct.fname = char(filename.Value);
end

%Read and display file description property
filedesclenIn = 0;
%Get the length of the 'DDC_FILE_DESCRIPTION' string property
[err,dummyVar,filedesclen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_DESCRIPTION,filedesclenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfiledesc=libpointer('stringPtr',blanks(filedesclen));
    [err,dummyVar,filedesc]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_DESCRIPTION,pfiledesc,filedesclen+1);
    setdatatype(filedesc,'int8Ptr',1,filedesclen);
    struct.desc = char(filedesc.Value);
end

%Read and display file title property
filetitlelenIn = 0;
%Get the length of the 'DDC_FILE_TITLE' string property
[err,dummyVar,filetitlelen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_TITLE,filetitlelenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfiletitle=libpointer('stringPtr',blanks(filetitlelen));
    [err,dummyVar,filetitle]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_TITLE,pfiletitle,filetitlelen+1);
    setdatatype(filetitle,'int8Ptr',1,filetitlelen);
    struct.title = char(filetitle.Value);
end

%Read and display file author property
fileauthlenIn = 0;
%Get the length of the 'DDC_FILE_AUTHOR' string property
[err,dummyVar,fileauthlen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_AUTHOR,fileauthlenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfileauth=libpointer('stringPtr',blanks(fileauthlen));
    [err,dummyVar,fileauth]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_AUTHOR,pfileauth,fileauthlen+1);
    setdatatype(fileauth,'int8Ptr',1,fileauthlen);
    struct.author = char(fileauth.Value);
end

% Get all the properties.
numFPropsIn = 0;
[err,numFProps] = calllib(libname, 'DDC_GetNumFileProperties', file, numFPropsIn);
if (err==0 && numFProps >= 1)
    struct.props = cell(numFProps, 2);
    for i=1:numFProps
        % Get file property name
        lenIn = 0;
        [err, len] = calllib(libname, 'DDC_GetFilePropertyNameLengthFromIndex', file, i-1, lenIn);
        
        if(err ~= 0)
            continue;
        end
        
        fpname = libpointer('stringPtr', blanks(len));
        [err, struct.props{i, 1}] = calllib(libname, 'DDC_GetFilePropertyNameFromIndex', file, i-1, fpname, len+1);
        
        if(err ~= 0)
            continue;
        end       
       
        % Now get its type.
        typeIn = 0;
        [err, ~, type] = calllib(libname, 'DDC_GetFilePropertyType', file, struct.props{i, 1}, typeIn);
               
       % If it's a double.
       typePtr = 'voidPtr';
       funcName = 'DDC_GetFileProperty';
       
       if strcmp(type,'DDC_Double')
          typePtr = 'doublePtr';
          funcName = 'DDC_GetFilePropertyDouble';
       elseif strcmp(type, 'DDC_Float')
          typePtr = 'singlePtr';
          funcName = 'DDC_GetFilePropertyFloat';
       elseif  strcmp(type, 'DDC_Int16')
           typePtr = 'int16Ptr';
           funcName = 'DDC_GetFilePropertyInt16';
       elseif strcmp(type, 'DDC_Int32')
           typePtr = 'int32Ptr';
           funcName = 'DDC_GetFilePropertyInt32';
       elseif strcmp(type, 'DDC_UInt8')
           typePtr = 'uint8Ptr';
           funcName = 'DDC_GetFilePropertyUInt8';
       elseif strcmp(type, 'DDC_String')
           % Need the string length
           lenIn = 0;
           [err, ~, len] = calllib(libname, 'DDC_GetFileStringPropertyLength', file, struct.props{i, 1}, lenIn);
           
           if(err ~= 0)
               continue;
           end
           
           % Get the value
           valIn = libpointer('stringPtr', blanks(len));
           val = libpointer('stringPtr', blanks(len));
           [err, ~, val] = calllib(libname, 'DDC_GetFileProperty', file, struct.props{i, 1}, valIn, len+1);
           setdatatype(val,'int8Ptr',1,len);
           if(err ~= 0)
               continue;
           end
           
           % Put the value in the struct
           struct.props{i, 2} = char(val.Value);    % 2 -> Value
           continue;  %We're done.
       elseif strcmp(type, 'DDC_Timestamp')
            % Not sure what to do with this yet.
           continue;          
       else
           continue;
       end
       
       val = libpointer('voidPtr', 0);
       valIn = libpointer('voidPtr', 0);
       [err, ~, val] = calllib(libname, funcName, file, fpname, valIn);
       setdatatype(val, typePtr, 1, 1);
       
       if(err == 0)
           struct.props{i, 2} = (val.Value);    %2 -> Value
       end
    end
end

%Read and display file timestamp property
yearIn = 0;
monthIn = 0;
dayIn = 0;
hourIn = 0;
minuteIn = 0;
secondIn = 0;
msecondIn = 0;
wkdayIn = 0;
[err,dummyVar,year,month,day,hour,minute,second,msecond,wkday]=calllib(libname,'DDC_GetFilePropertyTimestampComponents',file,DDC_FILE_DATETIME,yearIn,monthIn,dayIn,hourIn,minuteIn,secondIn,msecondIn,wkdayIn);
if err==0 %Only proceed if the property is found
    struct.timestamp = ['File Timestamp: ' num2str(month) '/' num2str(day) '/' num2str(year) ', ' num2str(hour) ':' num2str(minute) ':' num2str(second) ':' num2str(msecond)];
end

%Get channel groups
%Get the number of channel groups
numgrpsIn = 0;
[err,numgrps]=calllib(libname,'DDC_GetNumChannelGroups',file,numgrpsIn);
%Get channel groups only if the number of channel groups is greater than zero
if numgrps>0
	%Initialize an array to hold the desired number of groups
    pgrps=libpointer('int64Ptr',zeros(1,numgrps));
    [err,grps]=calllib(libname,'DDC_GetChannelGroups',file,pgrps,numgrps);
end

k = 1;
for i=1:numgrps %For each channel group
    %Get channel group name property
    grpnamelenIn = 0;
    [err,dummyVar,grpnamelen]=calllib(libname,'DDC_GetChannelGroupStringPropertyLength',grps(i),DDC_CHANNELGROUP_NAME,grpnamelenIn);
    if err==0 %Only proceed if the property is found
		%Initialize a string to the length of the property value
        pgrpname=libpointer('stringPtr',blanks(grpnamelen));
        [err,dummyVar,grpname]=calllib(libname,'DDC_GetChannelGroupProperty',grps(i),DDC_CHANNELGROUP_NAME,pgrpname,grpnamelen+1);
        setdatatype(grpname,'int8Ptr',1,grpnamelen);
    else
        grpname=libpointer('stringPtr','');
    end
    
    % Only want to read out one of these -> give them names if it's one of
    % them.
    name = char(grpname.Value);
      
    %Get channel group description property % (Unused)
    grpdesclenIn = 0;
    [err,dummyVar,grpdesclen]=calllib(libname,'DDC_GetChannelGroupStringPropertyLength',grps(i),DDC_CHANNELGROUP_DESCRIPTION,grpdesclenIn);
    if err==0 %Only proceed if the property is found
		%Initialize a string to the length of the property value
        pgrpdesc=libpointer('stringPtr',blanks(grpdesclen));
        [err,dummyVar,grpdesc]=calllib(libname,'DDC_GetChannelGroupProperty',grps(i),DDC_CHANNELGROUP_DESCRIPTION,pgrpdesc,grpdesclen+1);
    end
    
    %Get channels
    numchansIn = 0;
    %Get the number of channels in this channel group
    [err,numchans]=calllib(libname,'DDC_GetNumChannels',grps(i),numchansIn);
    %Get channels only if the number of channels is greater than zero
    if numchans>0
		%Initialize an array to hold the desired number of channels
        pchans=libpointer('int64Ptr',zeros(1,numchans));
        [err,chans]=calllib(libname,'DDC_GetChannels',grps(i),pchans,numchans);
    end
    
    cgstruct.channames=cell(1,numchans);
   
    % Get the channel group properties.
    npropin = 0;
    [err,nprop] = calllib(libname,'DDC_GetNumChannelGroupProperties', grps(i), npropin);
    if err == 0 % Only proceed if there are no errors.
        cgstruct.props = cell(nprop, 2);    % Dim 1-> Name, Dim2->Value
        for j = 1:nprop
            cpnamelenIn = 0;
            % Get the property name length.
            [err, cpnamelen] = calllib(libname, 'DDC_GetChannelGroupPropertyNameLengthFromIndex', grps(i), j-1, cpnamelenIn);
            if(err ~= 0 || cpnamelen < 1)
                continue;
            end
            
            cpnameIn = libpointer('stringPtr', blanks(cpnamelen));
            [err, cgstruct.props{j, 1}] = calllib(libname, 'DDC_GetChannelGroupPropertyNameFromIndex', grps(i), j-1, cpnameIn, cpnamelen+1);
            
            if(err ~= 0)
                continue;
            end
            
            % Now get the property's value
            
           % Start with its type.
           typeIn = 0;
           [err, ~, type] = calllib(libname, 'DDC_GetChannelGroupPropertyType', grps(i), cgstruct.props{j, 1}, typeIn);

           % If it's a double.
           typePtr = 'voidPtr';
           funcName = 'DDC_GetChannelGroupProperty';
           len = 1;
           
           if strcmp(type,'DDC_Double')
              typePtr = 'doublePtr';
              funcName = 'DDC_GetChannelGroupPropertyDouble';
              len = 8;
           elseif strcmp(type, 'DDC_Float')
              typePtr = 'singlePtr';
              funcName = 'DDC_GetChannelGroupPropertyFloat';
              len = 4;
           elseif  strcmp(type, 'DDC_Int16')
               typePtr = 'int16Ptr';
               funcName = 'DDC_GetChannelGroupPropertyInt16';
               len = 2;
           elseif strcmp(type, 'DDC_Int32')
               typePtr = 'int32Ptr';
               funcName = 'DDC_GetChannelGroupPropertyInt32';
               len = 4;
           elseif strcmp(type, 'DDC_UInt8')
               typePtr = 'uint8Ptr';
               funcName = 'DDC_GetChannelGroupPropertyUInt8';
               len = 1;
           elseif strcmp(type, 'DDC_String')
               % Need the string length
               lenIn = 0;
               [err, ~, len] = calllib(libname, 'DDC_GetChannelGroupStringPropertyLength', grps(i), cgstruct.props{j, 1}, lenIn);

               if(err ~= 0)
                   continue;
               end

               % Get the value
               valIn = libpointer('stringPtr', blanks(len));
               val = libpointer('stringPtr', blanks(len));
               [err, ~, val] = calllib(libname, 'DDC_GetChannelGroupPropertyString', grps(i), cgstruct.props{j, 1}, valIn, len+1);

               if(err ~= 0)
                   continue;
               end

               % Put the value in the struct
               cgstruct.props{j, 2} = val;    % 2 -> Value
               continue;  %We're done.
           elseif strcmp(type, 'DDC_Timestamp')
                % Not sure what to do with this yet.
               continue;          
           else
               continue;
           end

           valIn = libpointer(typePtr, 0);
           [err, ~, val] = calllib(libname, funcName, grps(i), cgstruct.props{j, 1}, valIn);
         
           if(err == 0)
               cgstruct.props{j, 2} = val;    %2 -> Value
           end

        end
    end
    
    for j=1:numchans %For each channel in the channel group
        %Get channel name property
        channamelenIn = 0;
        [err,dummyVar,channamelen]=calllib(libname,'DDC_GetChannelStringPropertyLength',chans(j),DDC_CHANNEL_NAME,channamelenIn);
        if err==0 %Only proceed if the property is found
			%Initialize a string to the length of the property value
            pchanname=libpointer('stringPtr',blanks(channamelen));
             [err,dummyVar,channame]=calllib(libname,'DDC_GetChannelProperty',chans(j),DDC_CHANNEL_NAME,pchanname,channamelen+1);
            setdatatype(channame,'int8Ptr',1,channamelen);
            cgstruct.channames{j}=char(channame.Value);
        else
            cgstruct.channames{j}='';
        end
        
        %Get channel data type
        typeIn = 0;
        [err,type]=calllib(libname,'DDC_GetDataType',chans(j),typeIn);
        
        %Get channel values for any given data type.
        
        types = {'DDC_Double', 'DDC_Float', 'DDC_Int16', 'DDC_Int32', 'DDC_UInt8'};
        typePtrs = {'doublePtr', 'singlePtr', 'int16Ptr', 'int32Ptr', 'uint8Ptr'};
        ind = find(strcmp(type, types));
        
        if ~isempty(ind)
            typePtr = typePtrs{ind};
            
            numvalsIn = 0;
            [err,numvals]=calllib(libname,'DDC_GetNumDataValues',chans(j),numvalsIn);
			%Initialize an array to hold the desired number of values
            
            pvals=libpointer(typePtr,zeros(1,numvals));
            [err,vals]=calllib(libname,'DDC_GetDataValues',chans(j),0,numvals,pvals);
            setdatatype(vals, typePtr, numvals)
            
            %Add channel values to a matrix. The comment, #ok<AGROW>, at
            %the end of the line prevents warnings about the matrix needing 
            %to allocate more memory for the added values.
            cgstruct.chanvals{j, :}=(vals.Value); %#ok<AGROW>
        else
            cgstruct.chanvals{j, :} = 0;  %#ok<AGROW> % Placeholder so that the channames and chanvals match up.
        end
            
    end    
   
    if(strcmp(name,MCTD_MAINDATA))
        cgstruct.name = MCTD_MAINDATA;
        struct.mdata = cgstruct;
    elseif(strcmp(name, MCTD_AVGDATA))
        cgstruct.name = MCTD_AVGDATA;
        struct.adata = cgstruct;
    elseif(strcmp(name, MCTD_PGNAME))
        cgstruct.name = MCTD_PGNAME;
        struct.prog = cgstruct;
    else
        struct.cgs{k} = cgstruct;
        k = k+1;
    end  
    
    clear cgstruct;
end

%Close file
err = calllib(libname,'DDC_CloseFile',file);

%Unload nilibddc.dll
unloadlibrary(libname);

% Now we want to go through the properties of struct.prog and turn it into
% something easily readable. We'll do the same for struct.mdata and
% struct.adata as well. b

% First find all the data we want and plug it into a struct called prog.

% Number of points
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PNP));
if isempty(ind)
   return; 
else
    prog.np = double(struct.prog.props{ind, 2});
    if(isempty(prog.np) || prog.np < 2)
        return;
    end
end

% Sample rate
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PSR));
if isempty(ind)
   return; 
else
    prog.sr = double(struct.prog.props{ind, 2});
    if(isempty(prog.sr) || prog.sr < 2)
        return;
    end
end

% Make a time vector
struct.t = (1:prog.np)/prog.sr;

% Number of transients
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PNT));
if isempty(ind)
   return; 
else
    prog.nt = double(struct.prog.props{ind, 2});
    if(isempty(prog.nt) || prog.nt < 1)
        return;
    end
end

% Number of dimensions
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PNDIMS));
if isempty(ind)
   return; 
else
    prog.nd = double(struct.prog.props{ind, 2});
    if(isempty(prog.nd))
        return;
    end
end

% Number of cycles
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PNCYCS));
if isempty(ind)
   return; 
else
    prog.ncyc = double(struct.prog.props{ind, 2});
    if(isempty(prog.ncyc))
        return;
    end
end

% Transient acquisition mode
% 0 -> ID First, then transients
% 1 -> All transients first, then ID
% 2 -> Phase cycles first, then IDS.
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PTMODE));
if isempty(ind)
    return;
else
    prog.ptmode = struct.prog.props{ind, 2};
    if(isempty(prog.ptmode))
        return;
    end
    
    if(prog.ptmode == 2)
        prog.ptmode = 1;
    end
end

% Number of instructions
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PNINSTRS));
if isempty(ind)
    return;
else
    prog.ninstr = double(struct.prog.props{ind, 2});
    if(isempty(prog.ninstr) || prog.ninstr < 1)
        return;
    end
end

% Number of unique instructions
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PNUINSTRS));
if isempty(ind)
    return;
else
    prog.nUniqueInstrs = double(struct.prog.props{ind, 2});
    if(isempty(prog.ninstr) || prog.nUniqueInstrs < prog.ninstr)
        return;
    end
end

% Whether or not it's varied
ind = find(strcmp(struct.prog.props(:, 1), MCTD_PVAR));
if isempty(ind)
    return;
else
    prog.varied = struct.prog.props{ind, 2};
    if(isempty(prog.varied))
        return;
    end
end

% These instructions are only for if it is varied
if prog.varied
    % Number of varied instructions
    ind = find(strcmp(struct.prog.props(:, 1), MCTD_PNVAR));
    if isempty(ind)
        return;
    else
        prog.nVar = double(struct.prog.props{ind, 2});
        if(isempty(prog.nVar))
            return;
        end
    end
    
    % Whether or not there's a skip condition
    ind = find(strcmp(struct.prog.props(:, 1), MCTD_PSKIP));
    if isempty(ind)
        return;
    else
        prog.skip = struct.prog.props{ind, 2};
        if(isempty(prog.skip))
            return;
        end
    end
    
    % Max number of steps
    ind = find(strcmp(struct.prog.props(:, 1), MCTD_PMAX_N_STEPS));
    if isempty(ind)
        return;
    else
        prog.maxnsteps = double(struct.prog.props{ind, 2});
        if(isempty(prog.maxnsteps))
            return;
        end
    end
    
    % Real number of steps (including skips)
    ind = find(strcmp(struct.prog.props(:, 1), MCTD_PREAL_N_STEPS));
    if isempty(ind)
        return;
    else
        prog.realnsteps = double(struct.prog.props{ind, 2});
        if(isempty(prog.realnsteps))
            return;
        end
    end
    
    % Grab maxsteps
    ind = find(strcmp(struct.prog.channames, MCTD_PMAXSTEPS));
    if isempty(ind)
        return;
    else
        prog.maxsteps = [struct.prog.chanvals{ind, 1}];
        if(isempty(prog.realnsteps))
            return;
        end
    end
    
    % Grab all the instructions
    prog.instrs{:, 1} = -1*ones(1, prog.nUniqueInstrs);
    prog.instrs{:, 2} = zeros(1, prog.nUniqueInstrs);
    prog.instrs{:, 4} = -1*ones(1, prog.nUniqueInstrs);

    % Flags
    ind = find(strcmp(struct.prog.channames, MCTD_PROGFLAG));

    if ~isempty(ind)
        vals = double([struct.prog.chanvals{ind, 1}]);
        for i = 1:length(vals)
            prog.instrs{i, 1} = vals(i);
        end
    end
    
    %Instr_delay
    ind = find(strcmp(struct.prog.channames, MCTD_PROGTIME));
    if ~isempty(ind)
        vals = double([struct.prog.chanvals{ind, 1}]);
        for i = 1:length(vals)
            prog.instrs{i, 2} = vals(i);
        end
    end
    
    % Instruction
    ind = find(strcmp(struct.prog.channames, MCTD_PROGINSTR));
    if ~isempty(ind)
        for i = 1:length(struct.prog.chanvals{ind, 1})
            ins = 'NA';
            vals = [struct.prog.chanvals{ind, 1}];
            switch vals(i)
                case 0
                    ins = 'CONTINUE';
                case 1
                    ins = 'STOP';
                case 2
                    ins = 'LOOP';
                case 3
                    ins = 'END_LOOP';
                case 4
                    ins = 'JSR';
                case 5
                    ins = 'RTS';
                case 6
                    ins = 'BRANCH';
                case 7
                    ins = 'LONG_DELAY';
                case 8
                    ins = 'WAIT';
            end
            prog.instrs{i, 3} = ins;
        end
    else
        for i = 1:prog.nUniqueInstrs
            prog.instrs{i, 3} = 'Not Found';
        end
    end
    
    % Instruction data
    ind = find(strcmp(struct.prog.channames, MCTD_PROGID));
    if ~isempty(ind)
       vals = double([struct.prog.chanvals{ind, 1}]);
       for i = 1:length(vals)
           prog.instrs{i, 4} = vals(i);
       end
    end
end

%%% TODO: Finish this to actually pull out the program %%%

% Now do the main data.

% Number of channels
ind = find(strcmp(struct.mdata.props(:, 1), MCTD_NCHANS));
if isempty(ind)
   return; 
else
    prog.nchans = double(struct.mdata.props{ind, 2});
    if(isempty(prog.nchans) || prog.nchans < 1)
        return;
    end
end

% Current step (linear index)
ind = find(strcmp(struct.mdata.props(:, 1), MCTD_CSTEP));
if isempty(ind)
   return; 
else
    prog.cstep = double(struct.mdata.props{ind, 2});
    if(isempty(prog.cstep))
        return;
    end
end

% Current step string
ind = logical(strcmp(struct.mdata.props(:, 1), MCTD_CSTEPSTR));
prog.cstepstr = struct.mdata.props{ind, 2};
struct.prog = prog;

% We have basically what we need now.
% Maxsteps now needs to be of the form [num chans, num transients, dim steps]
if(prog.nd < 1)
    maxsteps = [prog.np, prog.nchans, prog.cstep+1];
    ams = [prog.np, prog.nchans];
else
    maxsteps = num2cell(prog.maxsteps);
    [ms{1:(length(prog.maxsteps)+1)}] = ind2sub([prog.nt, maxsteps{(prog.ncyc+1):end}], prog.cstep);
    maxsteps = [prog.np, prog.nchans, ms{:}]; % Ignore the cycle crap
    ams = [prog.np, prog.nchans, ms{2:end}];
end

% Pre-allocate our data variable
mdata = zeros(maxsteps);
if(isfield(struct, 'adata'))
    adata = zeros(ams);
end

if(prog.nd < 1)
    % Main data first.
    for i = 1:prog.nchans
        
        if(prog.cstep < prog.nt)
            nt = prog.cstep;
        else
            nt = prog.nt;
        end
        
        maindata = [struct.mdata.chanvals{i, 1}];
        for j = 1:nt
            mdata(:, j, i) = maindata(((j-1)*prog.np+1):((j)*prog.np));
        end
    
        % Then average data
        if(isfield(struct, 'adata'))
            avgdata = [struct.adata.chanvals{i, 1}];
            adata(:, i) = avgdata(:);
        end
    end
   
    struct.mdata = mdata; % Set it in the actual struct
    
    if(isfield(struct, 'adata'))
        struct.adata = adata;
    end        
else
    % Get the linear indexing size for ind2sub/sub2ind
    switch(prog.ptmode)
        case 0  % ID First
            lsize = [ms{2:end}, ms{1}];
        case 1  % Transients first
            if(ms{2} > 1)
                lsize = [prog.nt, ms{2:end}];
            else
                lsize = ms{:};
            end
        case 2
            lsize = [ms{:}];
    end
    
    lsizeavg = double([ms{2:end}]);
    lsize = double(lsize);
    avgallsteps = prod(lsizeavg);
    allsteps = prod(lsize);
    
    if(ms{1} < prog.nt) 
        allsteps = allsteps - prog.nt;
        allsteps = allsteps+ms{1};
    end

    % Main data first.
    for i = 1:prog.nchans
        maindata = [struct.mdata.chanvals{i, 1}];
        for j = 1:allsteps
            [inds{1:length(lsize)}] = ind2sub(lsize, j); % Get the indices.
            cind = cell(prog.nd+prog.ncyc+1, 1);
           
            % Transform the subscripts so that the output will be uniform.
            % Always of the form [nt, ndim...]
            switch(prog.ptmode)
                case 0
                    cind(1) = inds(prog.nd+1);
                    cind(2:end) = inds(1:prog.nd);
                case 1
                    cind(:) = inds(:);
%               case 2
%                     cind(1) = prod(inds(1:prog.ncyc));
%                     cind(1) = cind(1)*inds(end);
%                     cind(2:end) = inds((prog.ncyc+1):(end-1));
            end
            
            mdata(:, i, cind{:}) = maindata(((j-1)*prog.np + 1):(j*prog.np));
        end
        
        struct.mdata = mdata;
        clear inds;
        
        % Now the average data
        if(isfield(struct, 'adata'))
            avgdata = [struct.adata.chanvals{i, 1}];
            for j = 1:avgallsteps
                [inds{1:length(lsizeavg)}] = ind2sub(lsizeavg, j);
                    
                adata(:, i, inds{:}) = avgdata(((j-1)*prog.np+1):(j*prog.np));
            end
        struct.adata = adata;    
        end
        
 
    end
   
    struct.mdata = mdata; % Set it in the actual struct
    
    if(isfield(struct, 'adata'))
        struct.adata = adata;
    end    
end

if(isfield(struct, 'adata'))
    s = num2cell(size(struct.adata));
    
    % Move the channels to the very end
    if(length(s) > 2)
        struct.adata = permute(struct.adata, [1, 3:length(s), 2]);
        
        if(struct.prog.nchans  == 1)
            struct.adata = reshape(struct.adata, s{1}, s{3:end});
        end
    end
end

if(isfield(struct, 'mdata'))
    s = num2cell(size(struct.mdata));
    
    % Move the channels to the very end
    if(length(s) > 2)
        struct.mdata = permute(struct.mdata, [1, 3:length(s), 2]);
        
        if(struct.prog.nchans  == 1)
            struct.mdata = reshape(struct.mdata, s{1}, s{3:end});
        end
    end
end



