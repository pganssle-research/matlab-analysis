function dataStruct=importTDMS(TDMSfilename,SaveConvertedFile,DisplayContents)
% 
% This function imports TDMS files into MATLAB, preserving all properties
% and data values.
%
%   dataStruct=importTDMS(TDMSfilename,SaveConvertedFile,DisplayContents)
% 
% INPUTS:
%   TDMSfilename (optional) - a string or cell array of strings which
%       provides the full path to the file(s) to be imported. If omitted,
%       the user will be prompted to select files. If this input is
%       omitted, other inputs cannot be used.
%   SaveConvertedFile (optional) - true/false flag which indicates whether
%       or not a .mat file should be created for each TDMS file. If 'true',
%       a .mat file is created in the same folder as the TDMS file with the
%       same file name (and .mat extension). The default value is 'false' 
%       (don't save). If this input is omitted, the third input cannot be
%       used.
%   DisplayContnets (optional) - true/false flag which indicates whether
%       the property names and values should be printed to the screen. This
%       provides an overview of a file's structure. The default value is
%       'false' (don't display).
% 
% OUTPUTS:
%   dataStruct - a structure (or structure array if importing multiple TDMS
%       files) which contains the contents of the TDMS file(s)
% 
% STRUCTURE OF 'dataStruct':
%   - dataStruct contains a top-level field 'dataStruct.filenameTDMS' which
%       contains the full path to the TDMS file that was converted
%   - Top-level properties are contained in the field array
%       'dataStruct.property'. Each property contains the fields
%       'dataStruct.property(i).name' and 'dataStruct.property(i).value'
%       which provide the name and value for the ith property.
%   - Channel groups are contained in the field array 'dataStruct.group'.
%       Each channel group may contain properties (same format as at the
%       top-level). Channel groups may contain channels in the field array
%       'dataStruct.group(i).channel'. Each channel may contain properties
%       (again as above). 'dataStruct.group(i).channel(j).data' contains
%       the array of data values for the channel. Groups and channels also
%       contain the field .name which brings that property up to the
%       channel/group level.
% 
% 
% EXAMPLE CALL AND RESULTING STRUCTURE
% >> dataStruct=importTDMS({'C:\file1.tdms','C:\file2.tdms'},true,false);
% 
% This results in the creation of 'C:\file1.mat' and 'C:\file2.mat' which
% each contain a single variable: dataStruct(1) and dataStruct(2), 
% respectively. Nothing is printed to the screen (third input is false),
% but would be similar to the following detailed structure explanation:
% 
% dataStruct(1).property(1).name = 'name'
% dataStruct(1).property(1).value = 'my first TDMS file'
% dataStruct(1).property(2).name = 'datetime'
% dataStruct(1).property(2).value = '21 OCT 2010, 9:36:35.515625'
% dataStruct(1).property(3).name = 'config_file'
% dataStruct(1).property(3).value = 'C:\CALCS.ini'
% dataStruct(1).group(1).name = 'Streaming'
% dataStruct(1).group(1).property(1).name = 'name'
% dataStruct(1).group(1).property(1).value = 'Streaming'
% dataStruct(1).group(1).channel(1).name = 'Input 1'
% dataStruct(1).group(1).channel(1).property(1).name = 'name'
% dataStruct(1).group(1).channel(1).property(1).value = 'Input 1'
% dataStruct(1).group(1).channel(1).property(2).name = 'unit_string'
% dataStruct(1).group(1).channel(1).property(2).value = 'Volts'
% dataStruct(1).group(1).channel(1).property(3).name = 'wf_increment'
% dataStruct(1).group(1).channel(1).property(3).value = '0.001'
% dataStruct(1).group(1).channel(1).data = [3.2 2.5 3.4 ... ]
% dataStruct(1).group(1).channel(2).name = 'Input 2'
% dataStruct(1).group(1).channel(2).property(1).name = 'name'
% dataStruct(1).group(1).channel(2).property(1).value = 'Input 2'
% dataStruct(1).group(1).channel(2).property(2).name = 'unit_string'
% dataStruct(1).group(1).channel(2).property(2).value = 'Amps'
% dataStruct(1).group(1).channel(2).property(3).name = 'wf_increment'
% dataStruct(1).group(1).channel(2).property(3).value = '0.001'
% dataStruct(1).group(1).channel(2).data = [0.1 0.1 0.2 ... ]
% dataStruct(1).group(2).name = 'Events'
% dataStruct(1).group(2).property(1).name = 'name'
% dataStruct(1).group(2).property(1).value = 'Events'
% dataStruct(1).group(2).channel(1).name = 'Trigger'
% dataStruct(1).group(2).channel(1).property(1).name = 'name'
% dataStruct(1).group(2).channel(1).property(1).value = 'Trigger'
% dataStruct(1).group(2).channel(1).data = ['21 OCT 2010, 9:37:00.000025' ...]
% 
% dataStruct(2).property(1).name = 'name'
% dataStruct(2).property(1).value = 'my second TDMS file'
% dataStruct(2).group(1).name = 'Constants'
% dataStruct(2).group(1).property(1).name = 'name'
% dataStruct(2).group(1).property(1).value = 'Constants'
% dataStruct(2).group(1).channel(1).name = 'Untitled'
% dataStruct(2).group(1).channel(1).property(1).name = 'name'
% dataStruct(2).group(1).channel(1).property(1).value = 'Untitled'
% dataStruct(2).group(1).channel(1).data = ['A' 'B' 'C' ... ]
%

% Author: Michael Corbett
%         AFRL/RZPE
%         michael.corbett@wpafb.af.mil
%         937-255-8977
% Last modified: 14 DEC 2010
% Approved for Public Release on 13 JAN 2011 (Case # 88ABW-2011-0105)


%% Check if the paths to 'nilibddc.dll' and 'nilibddc_m.h' can be found. If not, prompt the user to browse to each of the files.
NI_TDM_DLL_Path=which('nilibddc.dll');
libname='nilibddc';
if isempty(NI_TDM_DLL_Path) || iscell(NI_TDM_DLL_Path),
    [dllfile,dllfolder]=uigetfile('*dll','Select nilibddc.dll');
    libname=strtok(dllfile,'.');
    NI_TDM_DLL_Path=fullfile(dllfolder,dllfile);
end
NI_TDM_H_Path=which('nilibddc_m.h');
if isempty(NI_TDM_H_Path) || iscell(NI_TDM_H_Path),
    [hfile,hfolder]=uigetfile('*h','Select nilibddc_m.h');
    NI_TDM_H_Path=fullfile(hfolder,hfile);
end

%% Load nilibddc.dll (be sure to call 'unloadlibrary(libname)' when finished using the library)
loadlibrary(NI_TDM_DLL_Path,NI_TDM_H_Path);

%% Process function inputs
switch nargin
    case 0
        %Prompt the user to browse to the path of the TDM or TDMS file to read
        [filepath,filefolder]=uigetfile({'*.tdm;*.tdms'},'Select a TDM or TDMS file','MultiSelect','on');
        if isempty(filepath),
            return;
        end
        if ~iscell(filepath),
            TDMSfiles=cellstr(fullfile(filefolder,filepath)); 
        else
            for i=1:numel(filepath),
                TDMSfiles{i}=fullfile(filefolder,filepath{i});
            end
        end
        SaveConvertedFile=false;
        DisplayContents=false;
        
    case 1
        if ~ischar(TDMSfilename) && ~iscell(TDMSfilename),
            e=errordlg('If specified, the function''s first input argument must be a character string or a cell of character strings','Invalid Input Argument');
            uiwait(e)
            return
        end
        TDMSfiles=cellstr(TDMSfilename);
        SaveConvertedFile=false;
        DisplayContents=false;
        
    case 2
        if ~ischar(TDMSfilename) && ~iscell(TDMSfilename),
            e=errordlg('If specified, the function''s first input argument must be a character string or a cell of character strings','Invalid Input Argument');
            uiwait(e)
            return
        end
        if ~islogical(SaveConvertedFile) && ~ismember(SaveConvertedFile,[0,1])
            e=errordlg('If specified, the function''s second input argument must be ''true'' or ''false''','Invalid Input Argument');
            uiwait(e)
            return
        end
        TDMSfiles=cellstr(TDMSfilename);
        DisplayContents=false;
        
    case 3
        if ~ischar(TDMSfilename) && ~iscell(TDMSfilename),
            e=errordlg('If specified, the function''s first input argument must be a character string or a cell of character strings','Invalid Input Argument');
            uiwait(e)
            return
        end
        if ~islogical(SaveConvertedFile) && ~ismember(SaveConvertedFile,[0,1])
            e=errordlg('If specified, the function''s second input argument must be ''true'' or ''false''','Invalid Input Argument');
            uiwait(e)
            return
        end
        if ~islogical(DisplayContents) && ~ismember(DisplayContents,[0,1])
            e=errordlg('If specified, the function''s third input argument must be ''true'' or ''false''','Invalid Input Argument');
            uiwait(e)
            return
        end
        TDMSfiles=cellstr(TDMSfilename);
        
    otherwise
        e=errordlg('The function requires no more than 3 input arguments','Too Many Input Arguments');
        uiwait(e)
        help ReadFile2
        return
        
end

%% Process each file
for currentfile=1:numel(TDMSfiles),
    Data_Path=TDMSfiles{currentfile};

    % Open the file (Always call 'DDC_CloseFile' when you are finished using a file)
    clear fileIn; fileIn = 0;
    [err,dummyVar,dummyVar,file]=calllib(libname,'DDC_OpenFileEx',Data_Path,'',1,fileIn);

    % Get all top level properties
    numOfFileProps            = 0;
    dummynumber               = 0;

    % Get number of top level properties
    [err,numOfFileProps]=calllib(libname,'DDC_GetNumFileProperties',file,numOfFileProps);

    % Initialize arrays to hold properties
    topLevelFilePropNameLength=cell(1,numOfFileProps);
    topLevelFilePropNames=cell(1,numOfFileProps);
    topLevelFilePropValueType=cell(1,numOfFileProps);
    topLevelFilePropValues=cell(1,numOfFileProps);

    % Collect property name lengths, names, value types, an values
    for i=0:numOfFileProps-1, 
        % Get property name length and add it to array
        nameLength=uint32(0);
        [err,nameLength]=calllib(libname,'DDC_GetFilePropertyNameLengthFromIndex',file,i,nameLength);
        topLevelFilePropNameLength{i+1}=nameLength;

        % Get property name and add it to array
        propName=blanks(nameLength+1);
        [err,propName]=calllib(libname,'DDC_GetFilePropertyNameFromIndex',file,i,propName,nameLength+1);
        topLevelFilePropNames{i+1}=propName;

        % Get property value type
        [err,dummyVar,valueType]=calllib(libname,'DDC_GetFilePropertyType',file,topLevelFilePropNames{i+1},dummynumber);
        topLevelFilePropValueType{i+1}=valueType;

        % Get property value and add to array
        switch topLevelFilePropValueType{i+1}
            case 'DDC_UInt8'
                [err,dummyVar,topLevelFilePropValues{i+1}]=calllib(libname,'DDC_GetFilePropertyUInt8',file,topLevelFilePropNames{i+1},uint8(dummynumber));

            case 'DDC_Int16'
                [err,dummyVar,topLevelFilePropValues{i+1}]=calllib(libname,'DDC_GetFilePropertyInt16',file,topLevelFilePropNames{i+1},int16(dummynumber));

            case 'DDC_Int32'
                [err,dummyVar,topLevelFilePropValues{i+1}]=calllib(libname,'DDC_GetFilePropertyInt32',file,topLevelFilePropNames{i+1},int32(dummynumber));

            case 'DDC_Float'
                [err,dummyVar,topLevelFilePropValues{i+1}]=calllib(libname,'DDC_GetFilePropertyFloat',file,topLevelFilePropNames{i+1},single(dummynumber));

            case 'DDC_Double'
                [err,dummyVar,topLevelFilePropValues{i+1}]=calllib(libname,'DDC_GetFilePropertyDouble',file,topLevelFilePropNames{i+1},double(dummynumber));

            case 'DDC_String'
                [err,dummyVar,nameLength]=calllib(libname,'DDC_GetFileStringPropertyLength',file,topLevelFilePropNames{i+1},nameLength);
                propVal=blanks(nameLength+1);
                [err,dummyVar,propVal]=calllib(libname,'DDC_GetFilePropertyString',file,topLevelFilePropNames{i+1},propVal,nameLength+1);
                topLevelFilePropValues{i+1}=propVal;

            case 'DDC_Timestamp'
                yearIn = 0; monthIn = 0; dayIn = 0; hourIn = 0; minuteIn = 0; secondIn = 0; msecondIn = 0; wkdayIn = 0;
                [err,dummyVar,year,month,day,hour,minute,second,msecond,wkday]=calllib(libname,'DDC_GetFilePropertyTimestampComponents',file,topLevelFilePropNames{i+1},yearIn,monthIn,dayIn,hourIn,minuteIn,secondIn,msecondIn,wkdayIn);
                topLevelFilePropValues{i+1} = processTimestamp(err,year,month,day,hour,minute,second,msecond);

            otherwise
                e=errordlg(['Unrecognized datatype for property: ' topLevelFilePropNames{i+1}],'Bad property datatype');
                uiwait(e)
                return
        end
    end % end top level property processing loop


    % Get channel groups
    dummynumber = 0;

    % Get the number of channel groups
    numgrps = 0;
    [err,numgrps]=calllib(libname,'DDC_GetNumChannelGroups',file,numgrps);

    % Get channel groups only if the number of channel groups is greater than zero
    if numgrps>0
        % Initialize an array to hold the desired number of groups
        groupHandles=int32(zeros(1,numgrps));
        [err,groupHandles]=calllib(libname,'DDC_GetChannelGroups',file,groupHandles,numgrps);
    end

    % Initialize arrays to hold properties
    groupPropNameLength=cell(1,numgrps);
    groupPropNames=cell(1,numgrps);
    groupPropValueType=cell(1,numgrps);
    groupPropValues=cell(1,numgrps);
    groupNumChannels=cell(1,numgrps);
    channelPropNameLength=cell(1,numgrps);
    channelPropNames=cell(1,numgrps);
    channelPropValueType=cell(1,numgrps);
    channelPropValues=cell(1,numgrps);
    channelDataType=cell(1,numgrps);
    channelDataNumVals=cell(1,numgrps);
    channelDataValues=cell(1,numgrps);

    for i=1:numgrps % For each channel group
        % Get number of group level properties
        numOfGroupProps = 0;
        [err,numOfGroupProps]=calllib(libname,'DDC_GetNumChannelGroupProperties',groupHandles(i),numOfGroupProps);

        % Initialize arrays to hold group properties
        groupPropNameLength{i}=cell(1,numOfGroupProps);
        groupPropNames{i}=cell(1,numOfGroupProps);
        groupPropValueType{i}=cell(1,numOfGroupProps);
        groupPropValues{i}=cell(1,numOfGroupProps);

        % Collect group property name lengths, names, value types, and values
        for j=0:numOfGroupProps-1, 
            % Get property name length and add it to array
            nameLength=uint32(0);
            [err,nameLength]=calllib(libname,'DDC_GetChannelGroupPropertyNameLengthFromIndex',groupHandles(i),j,nameLength);
            groupPropNameLength{i}{j+1}=nameLength;

            % Get property name and add it to array
            propName=blanks(nameLength+1);
            [err,propName]=calllib(libname,'DDC_GetChannelGroupPropertyNameFromIndex',groupHandles(i),j,propName,nameLength+1);
            groupPropNames{i}{j+1}=propName;

            % Get property value type
            [err,dummyVar,valueType]=calllib(libname,'DDC_GetChannelGroupPropertyType',groupHandles(i),groupPropNames{i}{j+1},dummynumber);
            groupPropValueType{i}{j+1}=valueType;

            % Get property value and add to array
            switch groupPropValueType{i}{j+1}
                case 'DDC_UInt8'
                    [err,dummyVar,groupPropValues{i}{j+1}]=calllib(libname,'DDC_GetChannelGroupPropertyUInt8',groupHandles(i),groupPropNames{i}{j+1},uint8(dummynumber));

                case 'DDC_Int16'
                    [err,dummyVar,groupPropValues{i}{j+1}]=calllib(libname,'DDC_GetChannelGroupPropertyInt16',groupHandles(i),groupPropNames{i}{j+1},int16(dummynumber));

                case 'DDC_Int32'
                    [err,dummyVar,groupPropValues{i}{j+1}]=calllib(libname,'DDC_GetChannelGroupPropertyInt32',groupHandles(i),groupPropNames{i}{j+1},int32(dummynumber));

                case 'DDC_Float'
                    [err,dummyVar,groupPropValues{i}{j+1}]=calllib(libname,'DDC_GetChannelGroupPropertyFloat',groupHandles(i),groupPropNames{i}{j+1},single(dummynumber));

                case 'DDC_Double'
                    [err,dummyVar,groupPropValues{i}{j+1}]=calllib(libname,'DDC_GetChannelGroupPropertyDouble',groupHandles(i),groupPropNames{i}{j+1},double(dummynumber));

                case 'DDC_String'
                    [err,dummyVar,nameLength]=calllib(libname,'DDC_GetChannelGroupStringPropertyLength',groupHandles(i),groupPropNames{i}{j+1},nameLength);
                    propVal=blanks(nameLength+1);
                    [err,dummyVar,propVal]=calllib(libname,'DDC_GetChannelGroupPropertyString',groupHandles(i),groupPropNames{i}{j+1},propVal,nameLength+1);
                    groupPropValues{i}{j+1}=propVal;

                case 'DDC_Timestamp'
                    yearIn = 0; monthIn = 0; dayIn = 0; hourIn = 0; minuteIn = 0; secondIn = 0; msecondIn = 0; wkdayIn = 0;
                    [err,dummyVar,year,month,day,hour,minute,second,msecond,wkday]=calllib(libname,'DDC_GetChannelGroupPropertyTimestampComponents',groupHandles(i),groupPropNames{i}{j+1},yearIn,monthIn,dayIn,hourIn,minuteIn,secondIn,msecondIn,wkdayIn);
                    groupPropValues{i}{j+1} = processTimestamp(err,year,month,day,hour,minute,second,msecond);

                otherwise
                    e=errordlg(['Unrecognized datatype for property: ' groupPropNames{i}{j+1}],'Bad property datatype');
                    uiwait(e)
                    return

            end
        end % end channel group property processing loop

        
        % Get channels
        numchans = 0;
        % Get the number of channels in this channel group
        [err,numchans]=calllib(libname,'DDC_GetNumChannels',groupHandles(i),numchans);
        groupNumChannels{i}=numchans;

        % Get channels only if the number of channels is greater than zero
        if numchans>0
            % Initialize an array to hold the desired number of channels
            chanHandles=int32(zeros(1,numchans));
            [err,chanHandles]=calllib(libname,'DDC_GetChannels',groupHandles(i),chanHandles,numchans);

            % Initialize arrays to hold properties and data values for each channel
            channelPropNameLength{i}=cell(1,numchans);
            channelPropNames{i}=cell(1,numchans);
            channelPropValueType{i}=cell(1,numchans);
            channelPropValues{i}=cell(1,numchans);
            channelDataType{i}=cell(1,numchans);
            channelDataNumVals{i}=cell(1,numchans);
            channelDataValues{i}=cell(1,numchans);

            for j=1:numchans % For each channel in the channel group
                % Get number of channel properties and number of data elements
                numOfChannelProps = 0;
                [err,numOfChannelProps]=calllib(libname,'DDC_GetNumChannelProperties',chanHandles(j),numOfChannelProps);

                % Initialize arrays to hold properties and data values for each channel
                channelPropNameLength{i}{j}=cell(1,numOfChannelProps);
                channelPropNames{i}{j}=cell(1,numOfChannelProps);
                channelPropValueType{i}{j}=cell(1,numOfChannelProps);
                channelPropValues{i}{j}=cell(1,numOfChannelProps);

                % Collect channel property name lengths, names, value types, and values
                for k=0:numOfChannelProps-1, 
                    % Get property name length and add it to array
                    nameLength=uint32(0);
                    [err,nameLength]=calllib(libname,'DDC_GetChannelPropertyNameLengthFromIndex',chanHandles(j),k,nameLength);
                    channelPropNameLength{i}{j}{k+1}=nameLength;

                    % Get property name and add it to array
                    propName=blanks(nameLength+1);
                    [err,propName]=calllib(libname,'DDC_GetChannelPropertyNameFromIndex',chanHandles(j),k,propName,nameLength+1);
                    channelPropNames{i}{j}{k+1}=propName;

                    % Get property value type
                    [err,dummyVar,valueType]=calllib(libname,'DDC_GetChannelPropertyType',chanHandles(j),channelPropNames{i}{j}{k+1},dummynumber);
                    channelPropValueType{i}{j}{k+1}=valueType;

                    % Get property value and add to array
                    switch channelPropValueType{i}{j}{k+1}
                        case 'DDC_UInt8'
                            [err,dummyVar,channelPropValues{i}{j}{k+1}]=calllib(libname,'DDC_GetChannelPropertyUInt8',chanHandles(j),channelPropNames{i}{j}{k+1},uint8(dummynumber));

                        case 'DDC_Int16'
                            [err,dummyVar,channelPropValues{i}{j}{k+1}]=calllib(libname,'DDC_GetChannelPropertyInt16',chanHandles(j),channelPropNames{i}{j}{k+1},int16(dummynumber));

                        case 'DDC_Int32'
                            [err,dummyVar,channelPropValues{i}{j}{k+1}]=calllib(libname,'DDC_GetChannelPropertyInt32',chanHandles(j),channelPropNames{i}{j}{k+1},int32(dummynumber));

                        case 'DDC_Float'
                            [err,dummyVar,channelPropValues{i}{j}{k+1}]=calllib(libname,'DDC_GetChannelPropertyFloat',chanHandles(j),channelPropNames{i}{j}{k+1},single(dummynumber));

                        case 'DDC_Double'
                            [err,dummyVar,channelPropValues{i}{j}{k+1}]=calllib(libname,'DDC_GetChannelPropertyDouble',chanHandles(j),channelPropNames{i}{j}{k+1},double(dummynumber));

                        case 'DDC_String'
                            [err,dummyVar,nameLength]=calllib(libname,'DDC_GetChannelStringPropertyLength',chanHandles(j),channelPropNames{i}{j}{k+1},nameLength);
                            propVal=blanks(nameLength+1);
                            [err,dummyVar,propVal]=calllib(libname,'DDC_GetChannelPropertyString',chanHandles(j),channelPropNames{i}{j}{k+1},propVal,nameLength+1);
                            channelPropValues{i}{j}{k+1}=propVal;

                        case 'DDC_Timestamp'
                            yearIn = 0; monthIn = 0; dayIn = 0; hourIn = 0; minuteIn = 0; secondIn = 0; msecondIn = 0; wkdayIn = 0;
                            [err,dummyVar,year,month,day,hour,minute,second,msecond,wkday]=calllib(libname,'DDC_GetChannelGroupPropertyTimestampComponents',chanHandles(j),channelPropNames{i}{j}{k+1},yearIn,monthIn,dayIn,hourIn,minuteIn,secondIn,msecondIn,wkdayIn);
                            channelPropValues{i}{j}{k+1} = processTimestamp(err,year,month,day,hour,minute,second,msecond);

                        otherwise
                            e=errordlg(['Unrecognized datatype for property: ' channelPropNames{i}{j}{k+1}],'Bad property datatype');
                            uiwait(e)
                            return

                    end
                end % end channel property processing loop

                % Get channel data type and number of values
                [err,valueType]=calllib(libname,'DDC_GetDataType',chanHandles(j),dummynumber);
                channelDataType{i}{j}=valueType;
                [err,numVals]=calllib(libname,'DDC_GetNumDataValues',chanHandles(j),dummynumber);
                channelDataNumVals{i}{j}=numVals;

                % Get channel data
                switch channelDataType{i}{j}
                    case 'DDC_UInt8'
                        dataValues=uint8(zeros(1,numVals));
                        [err,dataValues]=calllib(libname,'DDC_GetDataValuesUInt8',chanHandles(j), uint32(0), channelDataNumVals{i}{j},dataValues);
                        channelDataValues{i}{j}=dataValues;

                    case 'DDC_Int16'
                        dataValues=int16(zeros(1,numVals));
                        [err,dataValues]=calllib(libname,'DDC_GetDataValuesInt16',chanHandles(j), uint32(0), channelDataNumVals{i}{j},dataValues);
                        channelDataValues{i}{j}=dataValues;

                    case 'DDC_Int32'
                        dataValues=int32(zeros(1,numVals));
                        [err,dataValues]=calllib(libname,'DDC_GetDataValuesInt32',chanHandles(j), uint32(0), channelDataNumVals{i}{j},dataValues);
                        channelDataValues{i}{j}=dataValues;

                    case 'DDC_Float'
                        dataValues=single(zeros(1,numVals));
                        [err,dataValues]=calllib(libname,'DDC_GetDataValuesFloat',chanHandles(j), uint32(0), channelDataNumVals{i}{j},dataValues);
                        channelDataValues{i}{j}=dataValues;

                    case 'DDC_Double'
                        dataValues=double(zeros(1,numVals));
                        [err,dataValues]=calllib(libname,'DDC_GetDataValuesDouble',chanHandles(j), uint32(0), channelDataNumVals{i}{j},dataValues);
                        channelDataValues{i}{j}=dataValues;

                    case 'DDC_String'
                        % *** Warning *** If the dataValues variable is not initialized as a
                        % cell array with enough space for each element (the hardcoded
                        % 100 characters in the next line) this may cause errors. I couldn't
						% figure out a way to know the length of individual data elements of
						% type DDC_String
                        dataValues=cell(1,numVals); for tmp=1:double(numVals), dataValues{tmp}=blanks(100); end
                        [err,dataValues]=calllib(libname,'DDC_GetDataValuesString',chanHandles(j), uint32(0), channelDataNumVals{i}{j},dataValues);
                        channelDataValues{i}{j}=dataValues;

                    case 'DDC_Timestamp'
                        year = uint32(zeros(1,numVals)); month = uint32(zeros(1,numVals)); day = uint32(zeros(1,numVals)); hour = uint32(zeros(1,numVals));
                        minute = uint32(zeros(1,numVals)); second = uint32(zeros(1,numVals)); msecond = zeros(1,numVals); wkday = uint32(zeros(1,numVals));
                        [err,year,month,day,hour,minute,second,msecond,wkday]=calllib(libname,'DDC_GetDataValuesTimestampComponents',...
                            chanHandles(j),uint32(0),channelDataNumVals{i}{j},year,month,day,hour,minute,second,msecond,wkday);
                        channelDataValues{i}{j} = processTimestamp(err,year,month,day,hour,minute,second,msecond);

                    otherwise
                        e=errordlg('Unrecognized datatype for channel: ','Bad property datatype');
                        uiwait(e)
                        return

                end
            end % end loop of channels within a group
        end % end if group contains channels
    end % end channel group processing loop


    % Close file
    err = calllib(libname,'DDC_CloseFile',file);

    % Post process the data into a nice MATLAB structure and optionally save and display it
    dataStruct(currentfile)=postProcess(Data_Path, topLevelFilePropNames, topLevelFilePropValues, groupPropNames, groupPropValues, channelPropNames, channelPropValues, channelDataValues);
    if DisplayContents, printStruct(dataStruct(currentfile)); end
    if SaveConvertedFile, 
        dataStructure=dataStruct(currentfile);
        matPath=strrep(Data_Path,'.tdms','.mat');
        matPath=strrep(matPath,'.TDMS','.mat');
        save(matPath,'dataStructure');
    end

end % end loop that processes each TDMS file

%% Unload nilibddc.dll
unloadlibrary(libname);

% %%%%%%%%%%% Some leftover debugging code if needed...
% topLevelFilePropNames'
% topLevelFilePropValues'
% for i=1:numgrps,
%     for j=1:groupNumChannels{i},
%         channelPropNames{i}{j}
%         channelPropValues{i}{j}
%         % channelDataValues{i}{j}
%     end
% end
% %%%%%%%%%%% End of debugging code...
end % End of main importTDMS function




% Subfunction to assist in postprocessing the data from several variables into one structure
function retStruct=postProcess(Data_Path, topLevelFilePropNames, topLevelFilePropValues, groupPropNames, groupPropValues, channelPropNames, channelPropValues, channelDataValues)
    % Create a structure variable
    retStruct=struct;
    
    % Set the TDMS filename as a top level field in the MATLAB structure
    retStruct.filenameTDMS = Data_Path;
    
    % Add all top level properties to the MATLAB structure as index .name and .value pairs
    for i=1:numel(topLevelFilePropNames), % for all top level properties
        retStruct.property(i).name=topLevelFilePropNames{i};
        retStruct.property(i).value=topLevelFilePropValues{i};
    end
    
    % Work through each of the channel groups, adding the contents to the MATLAB structure
    for i=1:numel(groupPropNames), % for all groups
        % Local variable to store channel group name
        clear groupName;
        % Find channel group name
        for j=1:numel(groupPropNames{i}), 
            if strcmp(groupPropNames{i}{j},'name'),
                groupName=groupPropValues{i}{j};
                break; % break for loop as soon as 'name' is found
            end
        end
        % If channel group name was not found, use this message as the name ...
        if ~exist('groupName', 'var'), groupName='Channel group name not found'; end;
        % Set the channel group name as a group level field
        retStruct.group(i).name=groupName;
        
        % Work through all the group properties, adding them as indexed .name and .value pairs
        for j=1:numel(groupPropNames{i}), % for all group properties
            retStruct.group(i).property(j).name=groupPropNames{i}{j};
            retStruct.group(i).property(j).value=groupPropValues{i}{j};
        end
        
        % Work through each of the channels, adding the contents to the MATLAB structure
        for j=1:numel(channelPropNames{i}), % for all channels
            % local variable to store channel name
            clear channelName;
            % Find channel name
            for k=1:numel(channelPropNames{i}{j}), % find and set channel name
                if strcmp(channelPropNames{i}{j}{k},'name'),
                    channelName=channelPropValues{i}{j}{k};
                    break; % break for loop as soon as 'name' is found
                end
            end
            % If channel name was not found, use this message as the name...
            if ~exist('channelName', 'var'), channelName='Channel name not found'; end;
            % Set the channel name as a channel level field
            retStruct.group(i).channel(j).name=channelName;
            
            % Work through all the channel properties, adding them as indexed .name and .value pairs
            for k=1:numel(channelPropNames{i}{j}), % for all channel properties
                retStruct.group(i).channel(j).property(k).name=channelPropNames{i}{j}{k};
                retStruct.group(i).channel(j).property(k).value=channelPropValues{i}{j}{k};
            end
            
            % Add actual channel data to the MATLAB structure
            retStruct.group(i).channel(j).data=channelDataValues{i}{j};
        end % end channel processing loop
    end % end channel group processing loop
end % end postProcess subfunction



% Subfunction to print structure contents to screen
function printStruct(retStruct)
    format compact;
    disp(' ');
    disp(['Contents of:  ' retStruct.filenameTDMS])
    disp('Top level properties:')
    
    % Print top level properties
    for i=1:numel(retStruct.property),
        disp([num2str(i) '   ' retStruct.property(i).name ' = ' num2str(retStruct.property(i).value)])
    end
    
    % Work through channel groups
    for i=1:numel(retStruct.group),
        disp(' ');disp(['Group ' num2str(i) ' (' retStruct.group(i).name ') properties:'])
        
        % Print group properties
        for j=1:numel(retStruct.group(i).property),
            disp([num2str(j) '   ' retStruct.group(i).property(j).name ' = ' num2str(retStruct.group(i).property(j).value)])
        end
        
        % Work through channels within the group
        for j=1:numel(retStruct.group(i).channel),
            disp(['Group ' num2str(i) ' (' retStruct.group(i).name '), Channel ' num2str(j) ' (' retStruct.group(i).channel(j).name ') properties:'])
            
            % Print channel properties
            for k=1:numel(retStruct.group(i).channel(j).property),
                disp([num2str(k) '   ' retStruct.group(i).channel(j).property(k).name ' = ' num2str(retStruct.group(i).channel(j).property(k).value)])
            end
        end % end loop of channels within a group
    end % end loop of channel groups
end % end subfunction printStruct



% Subfunction to assist in processing timestamp data
% Note that this function was designed to accept arrays for each of the inputs to
%    facilitate processing channel values that are all timestamps
function retStr=processTimestamp(err,year,month,day,hour,minute,second,msecond)
    if err==0,
        % Assume that this function was given arrays to process
        retStr=cell(1,numel(month));
        for i=1:numel(month)
            switch month(i)
                case 1,  monthTxt = 'JAN';
                case 2,  monthTxt = 'FEB';
                case 3,  monthTxt = 'MAR';
                case 4,  monthTxt = 'APR';
                case 5,  monthTxt = 'MAY';
                case 6,  monthTxt = 'JUN';
                case 7,  monthTxt = 'JUL';
                case 8,  monthTxt = 'AUG';
                case 9,  monthTxt = 'SEP';
                case 10, monthTxt = 'OCT';
                case 11, monthTxt = 'NOV';
                case 12, monthTxt = 'DEC';
            end
            retStr{i} = [num2str(day(i)) ' ' monthTxt ' ' num2str(year(i)) ', ' num2str(hour(i)) ':' num2str(minute(i)) ':' num2str(double(second(i))+msecond(i)/1000,8)];
        end
        
        % If the inputs were for a single timestamp rather than an array, fix the return type
        if numel(month) == 1, retStr=retStr{1}; end
        
    else % if there was an error when the timestamp components were gathered return an empty string
        retStr = '';
    end
end % end subfunction processTimestamp