function [out, prog, np, sr] = readout(chan, pathname)
	% Matlab Script for opening and processing file data from the console program
	% Specify the channel you want read out to the function. Multiple channel support to be added later.
	%
	% Usage:
	% [out, prog, np, sr] = readout(chan);
       
    if nargin < 1
        fprintf('No channel selected, default to 0\n')
        chan = 0;
    end
    
    if chan >= 8 || chan < 0
       fprintf('Invalid channel number.\n')
       return
    end
    
    % If you haven't supplied a pathname, we need to prompt for one
    hist = {};
    if nargin < 2
        if ~exist('readout_history.mat', 'file')
            default_dir = pwd; % If the readout_history file is missing, 
        else
            load('readout_history.mat');
            if length(hist) < 1 %#ok<*NODEF>
                default_dir = pwd;
            else
                % Now search through the history file to find the most
                % recently used one that exists. This is useful if the same
                % history file is synced across multiple systems. We'll set
                % a variable keepatmost in the readout_history.mat file, so
                % that we can adjust how long the history we want to keep
                % is. Default is keep all., keepatmost == -1 also means
                % keep all.
                for j = length(hist):-1:1 
                    if exist(hist{j}, 'file')
                        default_dir = hist{j};
                        break; % Stop looking once you've found it
                    elseif j == 1
                        default_dir = pwd;
                    end
                end
                dupes = ismember(hist, default_dir); % List of the positions of duplicate entries
                dupes = dupes(1:end-1); % The most recent one is OK to stay, the others shouldn't even be there.
                
                hist(dupes) = []; % Delete the relevant entries
            end
        end
        
        % If the better function is available, use it.
        if exist('uigetdir2.m', 'file')
            pathname = uigetdir2(default_dir); % If no pathname has been passed, prompt for one.
        else
            pathname = uigetdir(default_dir);
        end
    end
    
    
    if pathname == 0
        return
    end
    
    a = strfind(pathname, filesep); % Get the locations of all the path separators
    
    if a(end) == length(pathname)
        a = a(1:end-1); % If the last one is a path separator, we don't care about that
    end
    
    folder_name = pathname(a(end)+1:end); % Folder name, needed for all the stuff later
    
    % Now update the history
    dir_name = pathname(1:a(end)-1); % Directory that the folder is in
    dupes = ismember(hist, dir_name);
    hist(dupes) = [];
    
    hist = [hist, {dir_name}]; % Add the directory now
    
    if exist('keepatmost', 'var') && keepatmost >= 1 && length(hist) > keepatmost
        del_num = length(hist)-keepatmost;
        hist = hist((1+del_num):end);  %#ok<NASGU> % Delete the oldest entries first
    end
    
    if exist('readout_history.mat', 'file')
        save('readout_history.mat', 'hist', '-append'); % Save the readout history.
    else
        save('readout_history.mat', 'hist');
    end
   
   % The pathname should end in a filesep, the folder_name should not.
    if folder_name(end) == filesep
        folder_name = folder_name(1:end-1);
    end
    
    if pathname(end) ~= filesep
        pathname = strcat(pathname, filesep);
    end
    
    
	% The filename is always the folder name + .txt.
	% This is the base filename, and contains information about the entire experiment.
	file_name = strcat(folder_name, '.txt');
	filename = strcat(pathname, file_name);

	newfile = sprintf('%s%sChan%01d.mat', pathname, folder_name, chan); 
    if exist(newfile, 'file')
        resp = questdlg('Readout has already been performed and a file with the output has been generated. Overwrite this file, or generate a new one?', 'File already found', 'Load from File', 'Overwrite', 'Load from File');
        if strcmp(resp, 'Load from File')
            load(newfile);
            return
        end
    end
    
    if ~exist(filename, 'file')
        fprintf('File does not exist\n')
        fprintf('Path: %s\n', filename)
        return
    end

	% Get information about the experiment.
	channels = 1;
	num_points = 0;
	num_transients = 0;
	sampling_rate = 0.0;
    num_dim = 0;

	% This is all readout   
    % Import the data. If this is 1D, it will be a structure and the
    % headers will be in textdata. Otherwise it's a cell of just headers
    im = importdata(filename, '\n');
    
    if isstruct(im)
        im = im.textdata; % Throw away extra data if they are there
    end
    
    % Get rid of the unnecessary whitespace
    dels = ismember(im, '');
    im(dels) = [];
    
    for n = 1:length(im)
        % Stop and scan through these if this is of the form text= number
        line = im{n};
        [scan, count] = sscanf(im{n}, '%s %d'); 
        if count >=1
            % Check if it's the line for channels
			[scan, c] = sscanf(line, 'Channels: %d');
			if c > 0
				if scan <= 0
					fprintf('Channel number invalid\n')
					return
				end
				channels = scan;
				continue
			end

			% Number of points
			[scan, c] = sscanf(line, 'NPoints= %d');
			if c > 0 
				if scan < 1
					fprintf('Number of points invalid\n')
					return
				end
				num_points = scan;
				continue
			end

			% Sampling rate
			[scan, c] = sscanf(line, 'SamplingRate= %f');
			if c > 0
				if scan <= 0.0
					fprintf('Sampling rate invalid\n')
					return
				end
				sampling_rate = scan;
				continue
			end

			% Number of dimensions in the experiment
			[scan, c] = sscanf(line, 'nDimensions= %d');
			if c > 0
				if scan <= 0 || scan > 8
					fprintf('Number of Dimensions invalid\n')
					return
				end

				num_dim = scan;
				continue
            end
        end
        
  		% How many transients in each experiment
		[scan, count] = sscanf(line, 'TransientsCompleted= %*d of %d');
		if count <= 0
			continue
		end

		nt = scan(1);
        if nt >= 1
			num_transients = nt;
        end    
    end

	% Break on any situation where these guys were not found.
	if num_dim == 0 
		fprintf('Number of dimensions invalid/Not Found\n')
		return
	end

	if num_points == 0
		fprintf('Number of points not found\n')
		return
	end

	if num_transients == 0
		fprintf('Number of transients not found\n')
		return
	end

	if sampling_rate == 0.0
		fprintf('Sampling rate not found\n')
		return
	end

	max_points = [];

	% Variable to pass the relevant information to the appropriate function
	exp_info = [num_dim num_points num_transients channels];

    if num_dim > 1
		g = 0;
		% If it's multi-dimensional, we want to read out the number of indirect points in each dimension
		% We want to do this here before we pass it to the next function, because we already have the file open
  
        a = find(ismember(im, '[Point]'));
        
        if length(a) ~= 1
            fprintf('Indirect point matrix missing.\n')
            return
        end

		max_points = zeros(1, num_dim-1);
		for n=1:num_dim-1
			line = im{a+n};
			[scan, count] = sscanf(line, 'IndirectDim %d - %d of %*d');
			if count < 1
				fprintf('Invalid number of indirect points!\n')
				return
			end

			if scan(1) ~= n
				fprintf('Confused n-dimensional attributes file.\n')
				return
			end

			max_points(n) = scan(2);
		end
	end

	np = num_points;
	sr = sampling_rate; % Sample rate is in Hz.
	
	% Send it to the appropriate function depending on the number of dimensions.
	if num_dim > 1
		[out, prog] = get_data_nd(pathname, chan, max_points, exp_info);
	else
		[out, prog] = get_data_1d(pathname, chan, exp_info);
	end
	
	save(newfile, 'out', 'prog', 'np', 'sr'); % Save when done
	
	

function [transients, prog] = get_data_nd(pathname, chan, max_points, exp_info)
	% Reads out the n-dimensional data from a given experiment into an np x nd x nt array
	% Where np = number of points, nd = number of dimensions
	% Specify the channel you want to read out as an argument to the function
	%
	% Filenames constructed as follows:
	% Multichannel, multi-dimensional: basepathname\C0\D01-0001\D02-0001\...\D0n-0001\FIDTransient0001.txt

	% Read out the experiment info
	transients = [];
    prog = [];
    nd = exp_info(1);
	np = exp_info(2);
	nt = exp_info(3);
	c = exp_info(4);
    
      
	% Create the base path name
	basepath = pathname;
	if basepath(end) ~= filesep
		strcat(basepath, filesep);
    end

   
	if c > 1
		tstr = sprintf('C%01d%s', chan, filesep);
		basepath = strcat(basepath, tstr);
	end


	% Create transient output vector
	% Each column is a transient. 
	dims = [np nt max_points];
	trans_vec = zeros(dims);
		
	% Dummy check - this is how many indirect points there are in the experiment.
	fail_safe = prod(max_points);
    
    steps = fail_safe*nt; % Number of steps in the calculation
    step = 0; % The step you're on.
    
	current_point = ones(1, nd-1);
	first = 1;
	i = 0;

	n_varied = 0;
	varied_instrs = [];
	program = [];
	done = 0;
    
    h = waitbar(0, 'Readout 0% complete', 'Name', '0% complete', 'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1);');
    
    c = onCleanup(@()delete(h));
    
    setappdata(h, 'canceling', 0);
    old_perc = 0;
    while i < fail_safe
		i = i+1; 
		% Construct a base path name out of each point
		strbuf = '';
        for j=1:nd-1
			a = sprintf('D%02d-%04d%s', j, current_point(j), filesep);
			strbuf = strcat(strbuf, a);
        end
        
		c_dim_path = strcat(basepath, strbuf);
		for t = 1:nt
			% Iterate through the transients
			current_fname = sprintf('%sFIDTransient%04d.txt', c_dim_path, t);
   
			if ~exist(current_fname, 'file')
				fprintf('Invalid selection: %s\n', current_fname)
				return
			end
           
            % Add a progress bar
            step = step+1;
            perc = floor((step/steps)*100); % Calculate the percentage to the nearest percent
            
            % Update it only if it needs updating.
            if getappdata(h, 'canceling')
                delete(h);
                return
            end
            
            if perc > old_perc
                waitbar(perc/100, h, sprintf('Readout %d%% Complete', perc));
                set(h, 'Name', sprintf('%d%% Complete', perc));
            end
            
			% In the first file, we need to read out how the program goes
			% This allows us to use these as points in the indirect dimensions.          
            if first == 1
                [data, varied_instrs] = read_data(current_fname, 1);
                
     			% Pre-allocate the program array, then stop performing this readout.
				% Rows are the varied instruction
				% Column 1 is the instruction number
				% Column 2 is the delay in nanoseconds
				% The remaining matrices are the dimensions.
                n_varied = size(varied_instrs);
                n_varied = n_varied(1);
				s = [n_varied 2 max_points];
				program = zeros(s);

				first = 0;
            else
                data = read_data(current_fname);
            end

			if length(data) ~= np
				fprintf('Number of points does not match! Should be %d, but it is %d\n', np, length(data))
				return
			end

			% Fill out the transient data
			
			for p = 1:np
				cp = [p t current_point];
				cp = num2cell(cp);
				trans_vec(cp{:}) = data(p);
			end

			% Fill out the relevant program point]
			
			for l = 1:n_varied
				% Recall that varied_instrs has the following form:
				% 1: Original instruction number (probably useless)
				% 2: Initial value
				% 3: Increment value
				% 4: Dimension

				% Which instruction
				pp = [l 1 current_point];
				pp = num2cell(pp);
				program(pp{:}) = varied_instrs(l, 1);
				
				% What's the delay?
				pp = [l 2 current_point];
				pp = num2cell(pp);
				
				% Needs to be the point-1, because the mulitplier starts at 0.
				c_dim_point = current_point(varied_instrs(l, 4))-1;
				program(pp{:}) = varied_instrs(l, 2) + c_dim_point*varied_instrs(l, 3);
			end
		end

		% Increment to the next point.
		% This is done by iterating through the dimensions until you find the lowest significant
		% dimension that is not at the maximum, incrementing that and setting all the other bits
		% to 1. This is similar to how you would count:
		% 1999 you would see 9 == max, 9 == max, 9==max, 1 < max, so set the 1 bit to 2, and the others to 0 = 2000

		for d = 1:nd-1
			if current_point(d) < max_points(d)
				current_point(d) = current_point(d)+1;
				for k = 1:d-1
					current_point(k) = 1;
				end
				break
			end
		end

		% We're almost done!
		% Since this is a 1-based index, we need to go through one more time.
        if done == 1
            break
        end     
        
        if current_point == max_points
            done = 1;
			continue
        end
		
    end
    delete(h);

	% I set the output vectors here, rather than setting them during the program.
	% This is just in case I want to add some more error checking in later.
	transients = trans_vec;
	prog = program;
	
	
	
function [transients, program] = get_data_1d(pathname, chan, exp_info)
	% Reads out the 1-dimensional data from a given experiment into an np x nt array
	% Where np = number of points, nt = number of transients
	% Specify the channel you want to read out as an argument to the function
	%
	% Filenames constructed as follows:
	% Multichannel: basepathname\C0\basefoldername-FIDTransient%04d.txt

	% Read out the experiment info
	nd = exp_info(1);
	np = exp_info(2);
	nt = exp_info(3);
	c = exp_info(4);

	% Create the base path name
	basepath = pathname;
	if basepath(end) == filesep
		basepath = basepath(1:end-1);
	end

	% Get the folder name
	% Start with the basepath, which we know ends in the file separator (get rid of that)
	% To be error-resistant, we could mess with this more, but we'll assume no one will feed this horrible paths
	% Since this is a program only I'll be using mostly anyway.
	a = strfind(basepath(1:end), filesep);
	folder_name = basepath(a(end)+1:end);

	basepath = strcat(basepath, filesep);
	
	if c > 1
		tstr = sprintf('C%01d%c', chan, filesep);
		basepath = strcat(basepath, tstr);
	end

	trans_vec = zeros(np, nt);

	for i=1:nt
		current_fname = sprintf('%s%s-FIDTransient%04d.txt', basepath, folder_name, i);

		if ~exist(current_fname, 'file')
			fprintf('Malformed path!\n')
			return
        end
		
	
        if i == 1
            [data program] = read_data(current_fname, 2);
        else
            [data] = read_data(current_fname);
        end
                
		if length(data) ~= np
			fprintf('Wrong number of points!\n')
			return
		end

		trans_vec(:, i) = data;
	end

	transients = trans_vec;
		
    
	function [output, program] = read_data(filename, get_prog)
	% Reads data formatted appropriately.
    %
    % If get_prog == 0, don't read out the program at all
    % If get_prog == 1, read out the varied instructions in the 2D program
    % If get_prog == 2, read out the specific 1D program

    
    d = importdata(filename, '\n');
    
    output = d.data; % Super easy! Must check that this is the right size on return, though
      
    % Here's the hard part - taking this stuff and getting the program out
    % of it.
    if nargin < 2 || get_prog == 0
        program = [];
        return
    elseif get_prog == 1
        % If get_prog == 1, we need the varied instructions matrix
        
        % These give you the indices where you can find the number of
        % instructions and each instruction.
        nvind = find(~cellfun('isempty', strfind(d.textdata, 'nVaried=')));
        v_instr_inds = find(~cellfun('isempty', strfind(d.textdata, 'VaryInstr')));
       
        % Get the number of instructions, check that it's right
        [scan, count] = sscanf(d.textdata{nvind(1)}, 'nVaried= %d');
        if count ~= 1
            fprintf('Number of varied instructions not found in file %s.\n', filename);
            return
        end
        
        nv = scan(1);
        
        if nv ~= length(v_instr_inds)
            fprintf('Number of varied instructions does not match file in %s.\n', filename);
            return
        end
            
        % Make a matrix to contain this information.
        % Rows are instructions, the columns are as follows:
        % 1: Original instruction number (probably useless)
        % 2: Initial value
        % 3: Increment value
        % 4: Dimension
        varied_instrs = zeros(nv, 4);
        
        for n=1:nv
            % The two non-captured values are final and final_val, which are not needed.
            [scan, count] = sscanf(d.textdata{v_instr_inds(n)}, 'VaryInstr %d %d %f %f %f %f %*f %*f %d');
            if count > 0
                if nv == 0
                    fprintf('Invalid pulse program\n')
                    return
                end 
                % Val * val_units gives you the number of nanoseconds.
                v_instr = scan(1)+1;
                instr = scan(2);
                init_val = scan(3);
                init_units = scan(4);
                inc = scan(5);
                inc_units = scan(6);
                dim = scan(7);

                varied_instrs(v_instr, :) = [instr init_val*init_units inc*inc_units, dim];
            end
        end  
        
        program = varied_instrs;
    elseif get_prog == 2
        % If get_prog == 2, then we want the actual program for this
        % specific file.
        
        % Get the number of instructions and the instruction indices
        % Getting it twice is really just a check to make sure it's a good
        % file.
        ninstr_ind = find(~cellfun('isempty', strfind(d.textdata, 'NInstructions=')));
        instr_inds = find(~cellfun('isempty', strfind(d.textdata, 'Instruction ')));
        
        [scan, count] = sscanf(d.textdata{ninstr_ind(1)}, 'NInstructions= %d');
        if count ~= 1
            fprintf('Number of instructions not found in file %s.\n', filename);
            return
        end
        
        ninstr = scan(1);
        
        if ninstr ~= length(instr_inds)
            fprintf('Number of instructions listed doesn''t match number of instructions found in file %s.\n', filename);
            return
        end
        
        % Now create a program which is a structure containing the
        % following fields:
        instr_explain =  ['ttlarray => ni x 24 array of the logical values of the TTLS\n' ...
            'instr => ni x 1 array, values of the instruction (Continue, Stop, ETC).\n' ...
            'instr_data => ni x 1 array, value of the instruction data\n' ...
            'delay => ni x 1 array, delay times in nanoseconds\n'...
            'scan => ni x 1 array, whether or not a scan occured at each instruction\n'...
            'instr_explain => This instruction string.\n'];
        
        instrs = struct('ttlarray', zeros(ninstr, 24),...
            'instr', zeros(ninstr, 1), 'instr_data', zeros(ninstr, 1), ...
            'delay', zeros(ninstr, 1), 'scan', zeros(ninstr, 1), ...
            'instr_explain', instr_explain);
        
        for n = 1:ninstr
            [scan, count] = sscanf(d.textdata{instr_inds(n)}, 'Instruction %*d %d %d %d %d %f %f');
            if count < 1
                fprintf('Invalid instruction at instr number %d in file %s.\n', n, filename);
                return
            end
            
            instrs.scan(n) = scan(1);
            instrs.instr(n) = scan(3);
            instrs.instr_data(n) = scan(4);
            instrs.delay(n) = scan(5)*scan(6); % With unit conversion
            
            ttls = scan(2); % Convert this to the proper array
            ttls = fliplr(dec2bin(ttls, 24));
            instrs.ttlarray(n, :) = str2num(ttls(:));        
        end
        
        program = instrs;
    end
    
