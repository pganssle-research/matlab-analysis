%This script reads in Varian FDF files (2D and 3D)
%Use: type ReadVarFDF at the matlab command prompt
%Currently 3D data is displayed through the middle of the 3D dimension
%rank1 matrix1 storage1 bits1 type1
%rank is the dimension of the file
%matrix is the number of data points in each directions

%storage (float, double, int), etc (string)-so far all are float32 type
%bits is an integer
%type (real, imag, abs)
%clear fid
%fid=fopen([pathname filename], 'r','ieee-le');
%written by Lana Kaiser, Varian, Dec. 2009


global pathname image_space

if exist('pathname')
    if ischar(pathname)
    [filename,pathname]=uigetfile('*.*','Select FDF file',pathname);
    else
    [filename,pathname]=uigetfile('*.*','Select FDF file');
    end
else

   [filename,pathname]=uigetfile('*.*','Select FDF file');  
end

if filename ~=0
    %Open the input file
    filename=[pathname filename];
    [fid,msg]=fopen(filename,'r');
       if fid<0
        % There is an error
        str='File Open Failed'
       else
        filename
       end
end



endian = 'b'; % Big Endian from Sun spectrometers

%Start reading the file
counter_line=0;
while 1
    line=fgetl(fid)
    %Get the important parameters: rank, matrix, format:
    stryes=findstr(' rank ',line); %determine 2D or 3D
    if ~isempty(stryes)
    
    %dependent on having only one space between float and rank, but in some
    %cases there are two spaces
         %if strmatch('float  rank = ', line) 
    rank1_index=isstrprop(line, 'digit');
    rank1 = str2num(line(find(rank1_index>0)));
    
    
    end
    
    stryes=findstr(' bits ',line); 
    if ~isempty(stryes)
    %if strmatch('float  bits = ', line)
       eqsign=findstr(line,'=');
       scolon=findstr(line,';');
      
       bits1 = str2num(strtrim(line(eqsign+1:scolon-1)));
       
    end

    stryes=findstr(' matrix[] ',line); 
    if ~isempty(stryes)
    
    %if strmatch('float  matrix[] = ', line) %size of the total matrix
        if rank1==2; %this is 2D case
        
        matrix1_index=isstrprop(line, 'digit');
        comma=findstr(line,',');
        totaldigits = (find(matrix1_index>0));
        dim(1) = str2num(line(totaldigits(1):comma-1));
        dim(2) = str2num(line(comma+2:totaldigits(end)));
         
        elseif rank1==3;
          
         matrix1_index=isstrprop(line, 'digit');
        comma=findstr(line,',');
        totaldigits = (find(matrix1_index>0));
        dim(1) = str2num(line(totaldigits(1):comma(1)-1));
        dim(2) = str2num(line(comma(1)+2:comma(2)-1));
        
        
        dim(3) = str2num(line(comma(2)+2:totaldigits(end)))
        end
    end
    
    
    
     if strmatch('int    bigendian', line)
        [token, rem] = strtok(line,'int   bigendian = ');
        endian = 'l'; % PC format (linux) 
        if str2num(token)==1;
            endian = 'b';
        end
     
     end
     
     counter_line=counter_line+1;
    
    if isspace(line)|(counter_line >100),break,end
    
end

dim
endian
 %endian = 'l'; % PC format (linux) 


%We are done with figuring out the file size info, now move on to read the
%binary portion
if rank1==2; %this is 2D case
  disp('reading 2D...')
status = fseek(fid, -dim(1)*dim(2)*bits1/8, 'eof');

image_space=fread(fid,[dim(1), dim(2)],'float32',endian);
%a=fread(fid);
figure;pcolor(image_space);shading interp;colormap(gray);


else
%if rank1==3; 
%this is 3D case
    disp('reading 3D...')
    status = fseek(fid, -dim(1)*dim(2)*dim(3)*bits1/8, 'eof');
    
    %image_space=fread(fid,[dim(1), dim(2), dim(3)],'float32',endian);
    image_space=fread(fid,[dim(1)*dim(2), dim(3)],'float32',endian);
%a=fread(fid);
    image_space=reshape(image_space,[dim(1) dim(2) dim(3)]);
    NS=dim(3)/2;
figure;pcolor(squeeze(image_space(:,:,NS)));shading interp;colormap(gray);
end
    
if status<0
    disp('Error reading file')
end
%position=ftell(fid);
%image_space=fread(fid,[dim(1) dim(2)],'float32',endian);
%a=fread(fid);
%figure;pcolor(image_space);shading interp;colormap(gray);
fclose(fid);

