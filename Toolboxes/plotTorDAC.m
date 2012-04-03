function plotTorDAC(data,processNonStreaming)
% Plots all streaming channels of a TorDAC data file
% 
% plotTorDAC(data,processNonStreaming)
% 
% Just pass this function a processed TorDAC TMDS file's data structure. It
% will plot each channel in the 'Streaming' group vs. time.
% 
% Optionally pass a true/false flag to determine whether to plot data from
% the 'Point Log' channel group and print data from the 'Events' group to
% the screen (defaults to false).

% Author: Michael Corbett
%         AFRL/RZPE
%         michael.corbett@wpafb.af.mil
%         937-255-8977
% Last modified: 14 DEC 2010
% Approved for Public Release on 13 JAN 2011 (Case # 88ABW-2011-0105)

switch nargin
    case 1
        if ~isstruct(data),
            e=errordlg('The function''s first input argument must be a structure','Invalid Input Argument');
            uiwait(e)
            return
        end
        processNonStreaming = false;
    case 2
        if ~isstruct(data),
            e=errordlg('The function''s first input argument must be a structure','Invalid Input Argument');
            uiwait(e)
            return
        end
        if ~islogical(processNonStreaming) && ~ismember(processNonStreaming,[0,1]),
            e=errordlg('If specified, the function''s second input argument must be ''true'' or ''false''','Invalid Input Argument');
            uiwait(e)
            return
        end
    otherwise
        disp('Too many inputs specified');
        help plotTorDAC
        return
end

for i=1:numel(data.group), % for each channel group
    if strcmpi(data.group(i).name, 'Streaming'), % for the Streaming channel group
        for j=1:numel(data.group(i).channel), % for each channel in the group
            % find the channel units and dt
            clear t channelUnits channelIncrement;
            for k=1:numel(data.group(i).channel(j).property),
                if strcmp(data.group(i).channel(j).property(k).name,'Units'),
                    channelUnits=data.group(i).channel(j).property(k).value;
                elseif strcmp(data.group(i).channel(j).property(k).name,'wf_increment'),
                    channelIncrement=data.group(i).channel(j).property(k).value;
                end
            end % end finding channel details
            % check for missing channel details
            if ~exist('channelUnits', 'var'), channelUnits='Channel units not found'; end;
            if ~exist('channelIncrement', 'var'),
                disp(['The channel ''' data.group(i).channel(j).name ''' does not have ''wf_increment'' specified and cannot be plotted']);
            end;

            if exist('channelIncrement', 'var'), % plotable channel
                % create a figure using 'nextfig' if the function is available
                try, nextfig(data.group(i).channel(j).name); catch, figure; end;

                % create time vector
                t=0:channelIncrement:(numel(data.group(i).channel(j).data)-1)*channelIncrement;

                % plot variable, create title, and add axes lables
                plot(t,data.group(i).channel(j).data);
                grid minor; box on; xlabel('Time [s]');
                ylabel(['[' channelUnits ']']);
                title(data.group(i).channel(j).name,'Interpreter','none');
            end % end channel plot
        end % end Streaming group
        
    elseif processNonStreaming && strcmpi(data.group(i).name, 'Point Log'), % for the Point Log group
        for j=1:numel(data.group(i).channel), % find time log
            if strcmpi(data.group(i).channel(j).name,'time'),
                timeIndex=j;
            end
        end % end finding time log
        for j=1:numel(data.group(i).channel), % for each channel in the group
            if j~=timeIndex,
                % find the channel units
                clear channelUnits;
                for k=1:numel(data.group(i).channel(j).property),
                    if strcmp(data.group(i).channel(j).property(k).name,'Units'),
                        channelUnits=data.group(i).channel(j).property(k).value;
                    end
                end
                if ~exist('channelUnits', 'var'), channelUnits='Channel units not found'; end;
                % end finding channel units
                
                %create a figure using 'nextfig' if the function is available
                try, nextfig(data.group(i).channel(j).name); catch, figure; end;

                % plot variable, create title, and add axes lables
                hold all;grid minor;box on;
                for w=1:numel(data.group(i).channel(j).data),
                    plot(w,data.group(i).channel(j).data(w),'.','MarkerSize',25);
                end
                legend(data.group(i).channel(timeIndex).data,'Location','EastOutside')
                box on; xlabel('Log Point (date/time stamp in legend)');
                ylabel(['[' channelUnits ']']);
                title(data.group(i).channel(j).name,'Interpreter','none');
            end
        end % end Point Log Group

    elseif processNonStreaming && strcmpi(data.group(i).name, 'Events'), % for the Events group
        disp('Events:');
        for j=1:numel(data.group(i).channel), % find time log
            if strcmpi(data.group(i).channel(j).name,'time'),
                timeIndex=j;
            end
        end% end finding time log
        for j=1:numel(data.group(i).channel), % for each channel in the group
            if j~=timeIndex,
                for k=1:numel(data.group(i).channel(timeIndex).data),
                    disp([char(data.group(i).channel(timeIndex).data(k)) '  ==>  ' char(data.group(i).channel(j).data(k))]);
                end
            end
        end % end Events Group
    else
        disp(['The group ''' data.group(i).name ''' will not be processed']);
    end % end group name checking
end % end channel group loop