function s = parse_MestreC(file)
% Parses a MestreC file from .mrc file. This is essentially an XML file,
% but because it contains a binary container, xmlread can have problems
% with it if there are any 0x0 characters, so we'll, annoyingly, have to
% parse the whole thing ourselves.
%
% Currently this only supports NMR files, since this is the only thing I'm
% looking at.
%
% Returns a struct with all the parameters we've read out from the file.
%
% Usage:
% s = parse_mestreC(file);

if(~exist('file', 'var') || ~exist(file, 'file'))
	file = get_path('parce_mrc.mat', '.mrc');
end

f = fopen(file);

% Find the binary data, once that's been removed, we can feel free to throw
% away all the null values we want. They are in a tag called 'Points'
v_found = 0;
ptag = '<Points>';
ptage = '</Points>';
vtag = '<Values>';
vtage = '</Values>';
while(~feof(f))
	l = fgetl(f); % Read 100 characters from a line.
	% Find the <Values> tag.
	if ~isempty(strfind(deblank(l), vtag))
			v_found = 1;
			break;
	end
end

if ~v_found
	fclose(f);
	error('Values not found!');
end

% Find the <Points> tag, or the </Values> tag, whichever comes first.
while(~feof(f))
	pl = ftell(f); % Location of the points.
	l = fgetl(f);
	if ~isempty(strfind(deblank(l), ptag))
		break;
	end
	
	if ~isempty(strfind(deblank(l), vtage))
		fclose(f);
		error('Points not found.');
	end
end

fclose(f);

% As far as I know, the points are all on one line, so parse what's between
% <Points></Points>
llen = length(l); % How many bytes is the whole line.
ps = strfind(l, ptag);
pe = strfind(l, ptage);

f = fopen(file, 'rb');
ftn = 'Temp_file.xml'; % Copy everything but the byte array into an file.
ft = fopen(ftn, 'wb+');

sp = fread(f, pl); % Read everything up to the points array.
zero = sprintf('0');
sp(sp == 0) = zero;
fwrite(ft, sp);


fseek(f, pl+llen, 'bof'); % After the byte array
sp = fread(f);
fclose(f);

sp(sp == 0) = zero;
fwrite(ft, sp);
fclose(ft);

% OK, now we have our byte array and we have a sanitized XML file, we need
% to read it into a valid struct.
s = xml2struct(ftn);
delete(ftn); % It was a temporary file, as sad as we are to see it go.

s = s.NMRabc.Spectrum; % We don't care about version info.

% Parse the byte array now - we need some of these values to do it right.
endianness = 'n';
if(strcmpi('little', s.Main.Values.Endian.Text))
	endianness = 'l';
end

% Array size.
format = 'float32';
if(strcmpi('float32', s.Main.Values.Format.Text))
	format = 'float32';
	bs  = 4;
end

format = [format '=>double']; % Get our values as doubles.

f = fopen(file, 'rb');

nb = str2double(s.Main.Values.TotalBytes.Text);
bl = pe-(ps+length(ptag));
fstart = pl+ps+length(ptag) + (bl-nb);

fseek(f, fstart, 'bof');
barray = fread(f, nb/bs, format, 0, endianness);

fclose(f);

np = 0;
for i = 1:length(s.Main.Window)
	np = np + str2double(s.Main.Window{i}.Points.Text);
end

n = (nb/bs)/np;

s.Main.Values.Points = barray;

if(strcmpi('yes', s.Main.Values.Interleaved.Text))
	b2 = barray(1:2:end) + 1i*barray(2:2:end);
	s.Main.Values.Points = b2;
end
