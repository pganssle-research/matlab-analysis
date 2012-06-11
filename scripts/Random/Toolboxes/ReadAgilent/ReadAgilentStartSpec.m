function[final_kspace,np,ntraces,nblocks]= ReadAgilentStartSpec(nspec)
%This file reads in Varian raw fid file 
%Output:final_kspace,np,ntraces,nblocks sw sfrq

%Example:
%[final_kspace,np,ntraces,nblocks]= ReadVarStart(3);

global pathname  np ntraces nblocks arraydim seqcon sw sfrq 

if ischar(pathname)
    [filename,pathname]=uigetfile('*.*','Select fid file for spectroscopy data',pathname);
else
[filename,pathname]=uigetfile('*.*','Select fid file for the spectroscopy data');
end

if filename ~=0
    %Open the input file
    filename=[pathname filename];
    [fid,msg]=fopen(filename,'rt');
       if fid<0
        % There is an error
        str='File Open Failed'
       else
        filename
       end
else
    error('File Open Stopped')
    final_kspace=0;
    np=0;
    ntraces=0;
    nblocks=0;
end

 [kspace,np,ntraces,nblocks] = ReadAgilent(filename);
 
%seqcon='nccnn' %this is for gems
%seqcon='ncsnn'; %this is for sems
if(ispc)
filenameProc=[pathname '\procpar']
else
 filenameProc=[pathname 'procpar'] %this is linux
end
   [fid,msg]=fopen(filenameProc,'rt');
    if fid<0
        
     %There is an error
     str=['File Open Failed'];
     errordlg(str,title,'modal');
     
     else
   
     sw = ReadProcpar('sw',filenameProc)
     sfrq = ReadProcpar('sfrq',filenameProc)
     disp(sprintf('Reading in Agilent Procpar files from %s',filename));
     arraydim=ReadProcpar('arraydim',filenameProc)
     seqcon = ReadProcpar('seqcon',filenameProc)
     seqfil = ReadProcpar('seqfil',filenameProc)
     ns = ReadProcpar('ns',filenameProc)
    end  
    
    final_kspace=kspace;

    
    if isequal(seqfil,'csi3d')
        petable = ReadProcpar('petable',filenameProc)
        pelist= ReadProcpar('pelist',filenameProc)
        disp('This requires petable in the same folder as the fid data, if pelist doesnt exist')
       %   if size(pelist,2)>1
       %   table=pelist.';
        %  else
        
        s=importdata([pathname petable]);
        
        %[pathname petable '1']
        %global table1 table2 table3
         table1=s.data(1,:);
         table2=s.data(2,:);
         table3=s.data(3,:);
         % table1=textread([pathname petable '1'],'%d','delimiter',' ','headerlines',1);
         % table2=textread([pathname petable '2'],'%d','delimiter',' ','headerlines',1);
         % table3=textread([pathname petable '3'],'%d','delimiter',' ','headerlines',1);
         % end
          etl=1;
 
        
         slicedim=max(table3)+1
         ydim=max(table2)
         xdim=max(table1)
         final_kspace=zeros(slicedim,ydim,xdim,np/2);
        
          sstart=1;
         
   %    for slicecounter=1:slicedim
       for tablecounter=1:nblocks
          
                
                 sliceval=table3(tablecounter)+1;
                 yval=table2(tablecounter);
                 xval=table1(tablecounter);
                 
                 final_kspace(sliceval,yval,xval,:)=kspace(tablecounter,:);
                 
       
       end
       disp('Done reading 3D CSI data ...')
    end   %end of csi3d
 
   
    
 % return
    
        %read additional table for data reconstruction
%    if seqcon(2)&&seqcon(3)~='n'
        if ~isequal(seqfil,'csi3d')
        disp('this is 2D csi data')
        
        %rearrange the data
        if ns==1 %This is single slice or fid
        final_kspace=reshape(final_kspace,[nblocks np/2 ntraces]);
        else %This is compressed multislice seqcon='nccsn'
        final_kspace=reshape(final_kspace,[nblocks np/2 ns ntraces/ns]); 
        final_kspace0=zeros(size(final_kspace,3), size(final_kspace,1),size(final_kspace,2),size(final_kspace,4));
        final_kspace1=zeros(size(final_kspace,3), size(final_kspace,1),size(final_kspace,4),size(final_kspace,2));
        %Take care of the slices, which is still dim(3):
        size(final_kspace)
        size(final_kspace0)
        
        %figure;plot(squeeze(real(final_kspace(1,:,2,8))))
        
        for slices=1:size(final_kspace,3)
                    final_kspace0(slices,:,:,:)=final_kspace(:,:,slices,:);
        end
        
   %test     figure;plot(squeeze(real(final_kspace0(1,:,2,8))))
                
                for ima=1:size(final_kspace0,4)
                    final_kspace1(:,:,ima,:)=final_kspace0(:,:,:,ima);
                end
                
               % figure;plot(squeeze(real(final_kspace1(1,:,2,8))))
                
              final_kspace=final_kspace1;
        clear final_kspace1 final_kspace0
        
        %figure;plot(squeeze(real(final_kspace(1,:,2,8))))
        
         end  %end of compressed multislice or 1 slice
        
       % z1=fftshift(fftn(final_kspace));
       
        if ndims(final_kspace)<4
        dim1=size(final_kspace);
        final_kspace1=zeros(dim1(1), dim1(3),dim1(2));
          for i=1:dim1(2)
            final_kspace1(:,:,i)=final_kspace(:,i,:);
          end;
        final_kspace=final_kspace1;
        clear final_kspace1;
        end
        size(final_kspace)
        
       %  figure;pcolor((squeeze(abs(z1(:,100,:)))));shading flat;colormap(gray)
        else %this is 3d csi or spec
          if ndims(final_kspace)>3
           figure;pcolor(abs(fftshift(fftn(squeeze(final_kspace(3,:,:,3))))));colormap(gray);shading interp
          else 
          figure;plot(real(final_kspace(nspec,:)));
    title('This is FID #1')
    figure;plot(real(fliplr(fftshift(fft(final_kspace(nspec,:))))));
    title('This is SPECTRUM #1')
          end
    end
 