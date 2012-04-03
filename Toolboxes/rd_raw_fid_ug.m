function [fid, ref_scan_loc] = rd_raw_fid_ug(dn)

% return the k-space data stripped of the header 
% dn = directory containing the fid file
% fid(x,y,k,n,t) = k-space point corresponding to readout step x, line y,
%          slice k, segment n and time t
% checked with multi-slice, non-segmented data set, okay, 1/12
% checked with multi-slice, segmented data set, okay, 1/12
% checked with reference scan location, okay, 1/17

ppPath = [dn strcat(filesep,'procpar')];
fidPath = [dn strcat(filesep,'fid')];
fp = fopen( fidPath, 'r', 'ieee-be');
fprintf(2, 'Reading: %s\n', fidPath);

mh = getMainHdr(fp);

np   = getPPV( 'np', ppPath);
pss  = getPPV( 'pss', ppPath);
nv   = getPPV( 'nv', ppPath);
nseg = getPPV( 'nseg', ppPath);
ref_scan_loc = getPPV( 'ref_pos', ppPath) + 1;


dx = np / 2;
num_slices = size(pss, 2);
%
dy = mh.ntraces;
%
dt = mh.nblocks;
fprintf( 1, 'Dims: %d x %d x %d\n', dx, dy, dt);

fid = zeros(dx,dy,dt) + i*ones(dx,dy,dt);

% fftScale = dx * dy;

mainHdrSize = 32;
blkHdrSize = 28;
offArr = zeros(1, mh.ntraces * mh.nblocks);
fileOff = mainHdrSize;
idx = 1;
for blkIdx=1:mh.nblocks
       fileOff = fileOff + blkHdrSize;
          for trcIdx=1:mh.ntraces
                    offArr(idx) = fileOff;
                            idx = idx + 1;
                                  fileOff = fileOff + mh.tbytes;
          end
end

if (mh.ebytes == 2) 
        precision = 'int16';
else
        precision = 'float';
end

realIdx = 1:2:np;
imagIdx = realIdx + 1;

wait_read = waitbar(0,'Raw data reading....');

idx = 1;
for it=1:mh.nblocks,
    %    fprintf( 2, '\nTime: %d ', it);
    %    fprintf( 2, '.');
    
    waitbar(it/mh.nblocks)   

    for iy=1:mh.ntraces,
          fprintf('Reading Block %d Trace %d.\n',it,iy);
        if (0 ~= fseek( fp, offArr(idx), 'bof'))
            error('Seek failed');
        end
        trc = fread( fp, np, precision);
        fid(:,iy,it) = (trc(realIdx) + 1i*trc(imagIdx));
        idx = idx + 1;
    end
end

close(wait_read);

fprintf( 2, '\nDone\n' );

if nv > 1
     %fid = reshape(fid,[dx floor(nv/nseg) num_slices nseg dt]);
end

fclose(fp);