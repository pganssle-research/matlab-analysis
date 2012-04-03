function out = mc_read_data(path)

if(~exist('path', 'var'))
    path = -1;
end

% Get the raw structure
[s, f] = mc_read_bin(path, 'mc_read_data_hist.mat'); %#ok

out = [];
if(isempty(f))
    return;
end

% Separately process the groups
MCD_DATAHEADER = '[Data Header]';
MCD_DISPHEADER = '[Display Header]';
MCD_DATAGROUP = '[DataGroup]';
MCD_PROGHEADER = '[PulseProgram]';

% Data header should come first - That'll be the main portion of the
% structure - so those are top-level values.
s1 = find_struct_by_name(f, MCD_DATAHEADER);
if(isempty(s1))
    return;
end

% Data Structure Names
MCD_FNAME = 'filename';
MCD_ENAME = 'ExperimentName';
MCD_ENUM = 'ExperimentNum';
MCD_HASH = 'HashCode';
MCD_NCHANS = 'NumChans';
MCD_TSTART = 'TimeStarted';
MCD_TDONE = 'TimeDone';
MCD_CIND = 'CurrentIndex';

out.FileName = [];
out.ExperimentName = [];
out.ExperimentNum = [];
out.hash = [];
out.tstart = [];
out.tdone = [];
out.nc = 0;
out.cind = -1;

sb = find_struct_by_name(s1.data, MCD_FNAME);
if(~isempty(sb))
    fname = deblank(sb.data');
    li = find(fname == '\', 1, 'last');
    if(isempty(li) || li == length(fname))
        out.FileName = fname;
    else
        out.FileName = fname((li+1):end);
    end
end

sb = find_struct_by_name(s1.data, MCD_ENAME);
if(~isempty(sb))
    out.ExperimentName = deblank(sb.data');
end

sb = find_struct_by_name(s1.data, MCD_ENUM);
if(~isempty(sb))
    out.ExperimentNum = sb.data;
end

sb = find_struct_by_name(s1.data, MCD_HASH);
if(~isempty(sb))
    out.hash = sb.data;
end

sb = find_struct_by_name(s1.data, MCD_TSTART);
if(~isempty(sb))
    out.tstart = deblank(sb.data');
end

sb = find_struct_by_name(s1.data, MCD_TDONE);
if(~isempty(sb))
   out.tdone = deblank(sb.data'); 
end

sb = find_struct_by_name(s1.data, MCD_NCHANS);
if(~isempty(sb))
   out.nc = sb.data; 
end

sb = find_struct_by_name(s1.data, MCD_CIND);
if(~isempty(sb))
    out.cind = sb.data;
end

% Display header is next - we can just do a direct dump
out.disp = [];
[~, loc] = find_struct_by_name(f, MCD_DISPHEADER);
if(isfield(s, loc))
    out.disp = eval(['s.' loc]);
end


function [s loc] = find_struct_by_name(in, name)
% Find a struct from its .name parameter.
s = [];
loc = [];
flist = fieldnames(in);

for i = 1:length(flist)
    b = eval(['in.' flist{i} ';']);
    
    if(isfield(b, 'name') && strcmp(b.name, name))
        s = b;
        loc = [flist{i}];
        break;
    end
    
    if(isstruct(b.data))
        [s l] = find_struct_by_name(b.data, name);
        if(~isempty(s))
            loc = [flist{i} '.' l];
            break;
        end
    end
end
