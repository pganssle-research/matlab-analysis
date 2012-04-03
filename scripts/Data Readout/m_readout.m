function [out, prog, np, sr] = m_readout(chan)
	% Matlab Script for opening and processing file data from the console program
	% Specify the channel you want read out to the function. Multiple channel support to be added later.

	default_file = 'C:\Documents and Settings\OmegaUser\My Documents\Omega Project\Programs\Magnetometer Controller\Experiments\\';
    default_file = 'C:\Documents and Settings\OmegaUser\My Documents\Omega Project\Programs\Magnetometer Controller\Experiments\T1 Measurements\Decane\110504Decane0000\\';
		
	%[filename, pathname] = uigetfile('*.txt;*.*', 'Select base file', default_file);

	%filename = strcat(pathname,filename);
    
    if nargin > 1
        fprintf('Incorrect number of arguments\n')
        return
    end
    
    if nargin < 1
        fprintf('No channel selected, default to 0\n')
        chan = 0;
    end
    
    if chan >= 8 || chan < 0
       fprintf('Invalid channel number.\n')
       return
    end
    
	pathname = uigetdir(default_file);
    a = strfind(pathname, filesep);
    folder_name = pathname(a(end)+1:end);
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

	if ~exist(filename, 'file')
		fprintf('File does not exist\n')
        fprintf('Path: %s\n', filename)
        return
    end

	fileid = fopen(filename, 'r');
    if fileid < 0
        fprintf('Invalid file\n')
        return
    end

	% Get information about the experiment.
	num_dimensions = 0;
	channels = 1;
	num_points = 0;
	num_transients = 0;
	sampling_rate = 0.0;

	% This is all readout
	while ~feof(fileid)
		line = fgetl(fileid);
		if strcmp(line, '[FID]') || strcmp(line, '[TransientData]')
			break
		end

		% Stop and scan through these if this is of the form text= number
		[scan, count] = sscanf(line, '%s %d');
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
		frewind(fileid);
		g = 0;


		% If it's multi-dimensional, we want to read out the number of indirect points in each dimension
		% We want to do this here before we pass it to the next function, because we already have the file open
		while ~feof(fileid)
			line = fgetl(fileid);
			if strcmp(line, '[Point]')
				g = 1;
				break
            end
        end

		if g == 0
			fprintf('Invalid number of indirect points!\n')
			return
		end

		max_points = zeros(1, num_dim-1);
		for i=1:num_dim-1
			line = fgetl(fileid);
			[scan, count] = sscanf(line, 'IndirectDim %d - %d of %*d');
			if count < 1
				fprintf('Invalid number of indirect points!\n')
				return
			end

			if scan(1) ~= i
				fprintf('Confused n-dimensional attributes file.\n')
				return
			end

			max_points(i) = scan(2);
		end
	end

	fclose(fileid);

	np = num_points;
	sr = sampling_rate; % Sample rate is in Hz.
	
	% Send it to the appropriate function depending on the number of dimensions.
	if num_dim > 1
		[out, prog] = get_data_nd(pathname, chan, max_points, exp_info);
	else
		[out] = get_data_1d(pathname, chan, exp_info);
		prog = []
	end
	
	newfile = strcat(pathname, folder_name);
	newfile = strcat(newfile, '-Array.mat');
	save(newfile, 'out', 'prog', 'np', 'sr');
	
	

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
		strcat(basepath, tstr);
	end


	% Create transient output vector
	% Each column is a transient. 
	dims = [np nt max_points];
	trans_vec = zeros(dims);
		
	% Dummy check - this is how many indirect points there are in the experiment.
	fail_safe = 1;
	for j = 1:length(max_points)
		fail_safe = fail_safe * max_points(j);
	end

	current_point = ones(nd-1);
	first = 1;
	i = 0;

	n_varied = 0;
	varied_instrs = [];
	program = [];
	done = 0;
    
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
				fprintf('Invalid selection\n')
				return
			end

			fileid = fopen(current_fname, 'r');

			if fileid < 0
				fprintf('Invalid selection\n')
				return
			end


			% In the first file, we need to read out how the program goes
			% This allows us to use these as points in the indirect dimensions.
			if first == 1
				while ~feof(fileid)
					line = fgetl(fileid);
					% Number of instructions total.
					[scan, count] = sscanf(line, 'nVaried= %d');
					if count > 0
						if scan(1) < 1
							fprintf('Invalid pulse program\n')
							return
						end
						n_varied = scan(1);

						% Make a matrix to contain this information.
						% Rows are instructions, the columns are as follows:
						% 1: Original instruction number (probably useless)
						% 2: Initial value
						% 3: Increment value
						% 4: Dimension
						varied_instrs = zeros(n_varied, 4);
						continue
					end

					% The two non-captured values are final and final_val, which are not needed.
					[scan, count] = sscanf(line, 'VaryInstr %d %d %f %f %f %f %*f %*f %d');
					if count > 0
						if n_varied == 0
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
						continue
					end				

					% Break when we get to the end, before we read out the transient data..
					[scan, count] = sscanf(line, '[TransientData]');
					if count > 0
						break
					end
				end

				if n_varied == 0
					fprintf('Invalid Pulse Program\n')
					return
				end

				% Pre-allocate the program array, then stop performing this readout.
				% Rows are the varied instruction
				% Column 1 is the instruction number
				% Column 2 is the delay in nanoseconds
				% The remaining matrices are the dimensions.
				s = [n_varied 2 max_points];
				program = zeros(s);

				first = 0;
			end

			% Pass to the data readout file
			data = read_data(fileid);

			fclose(fileid);

			if length(data) ~= np
				fprintf('Number of points does not match! Should be %d, but it is %d\n', np, length(data)
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
		
		fprintf('Point= %d\n', current_point(1))
		
	end

	% I set the output vectors here, rather than setting them during the program.
	% This is just in case I want to add some more error checking in later.
	transients = trans_vec;
	prog = program;
	
	
	
function [transients] = get_data_1d(pathname, num_dim, c)
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
	if basepath(end) ~= '\\'
		strcat(basepath, '\\');
	end

	% Get the folder name
	folder = strsplit(pathname, pathsep);
	folder_name = folder(len(folder));

	if c > 1
		tstr = sprintf('C%01d\\', chan);
		strcat(basepath, tstr);
	end

	trans_vec = zeros(np, nt);

	for i=1:nt
		current_fname = sprintf('%s-FIDTransient%04d.txt', basepath, i);

		if ~exist(current_fname, 'file')
			fprintf('Malformed URL!\n')
			return
		end
		
	
		fileid = fopen(current_fname);
		data = read_data(fileid);
		fclose(fileid);

		if length(data) ~= np
			fprintf('Wrong number of points!\n')
			return
		end

		trans_vec(:, i) = data;
	end

	transients = trans_vec;
		

	function [output] = read_data(fileid)
	% Reads data formatted appropriately.

	% Change this to wherever the default directory is where you store files.
	% TODO: Implement a version of this stored in a preferences file

	if fileid == -1
		fprintf('Failed to open file\n')
		return
	end

	frewind(fileid)

	while ~feof(fileid)
		l = fgetl(fileid);
		if strcmp(l, '[FID]') || strcmp(l, '[TransientData]')
			break
		end
	end

	A = [];

	while ~feof(fileid)
		l = fgetl(fileid);
		l = str2num(l);
		if ~isscalar(l)
			break
		end

		A = [A; l];
	end

	output = A;
