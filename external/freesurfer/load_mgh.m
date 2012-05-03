function [vol, M, mr_parms, volsz] = load_mgh(fname,slices,frames,headeronly)
% [vol, M, mr_parms, volsz] = load_mgh(fname,<slices>,<frames>,<headeronly>)
%
% fname - path of the mgh file
% 
% slices - list of one-based slice numbers to load. All
%   slices are loaded if slices is not specified, or
%   if slices is empty, or if slices(1) <= 0.
%
% frames - list of one-based frame numbers to load. All
%   frames are loaded if frames is not specified, or
%   if frames is empty, or if frames(1) <= 0.
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
% indices are 0-based. If the input has multiple frames,
% only the first frame is read.
%
% mr_parms = [tr flipangle te ti fov]
%
% volsz = size(vol). Helpful when using headeronly as vol is [].
%
% See also: save_mgh, vox2ras_0to1
%


%
% load_mgh.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision$
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

vol = [];
M = [];
mr_parms = [];
volsz = [];

if(nargin < 1 | nargin > 4)
  msg = 'USAGE: [vol M] = load_mgh(fname,<slices>,<frames>,<headeronly>)';
  fprintf('%s',msg);
  return;
end

% unzip if it is compressed 
if (strcmpi(fname((strlen(fname)-3):strlen(fname)), '.MGZ') | ...
		strcmpi(fname((strlen(fname)-3):strlen(fname)), '.GZ'))
  rand('state', sum(100*clock));
  gzipped =  round(rand(1)*10000000 + ...
		   sum(int16(fname))) + round(cputime);
  ind = findstr(fname, '.');
  new_fname = sprintf('/tmp/tmp%d.mgh', gzipped);
  if (ismac || strcmp(computer,'MAC') || strcmp(computer,'MACI') || strcmp(computer, 'MACI64'))
    unix(sprintf('gunzip -c %s > %s', fname, new_fname)) ;
  else
    unix(sprintf('zcat %s > %s', fname, new_fname)) ;
  end
  fname = new_fname ;
else
  gzipped = -1 ;
end


if(exist('slices')~=1) slices = []; end
if(isempty(slices)) slices = 0; end
if(slices(1) <= 0) slices = 0; end

if(exist('frames')~=1) frames = []; end
if(isempty(frames)) frames = 0; end
if(frames(1) <= 0) frames = 0; end

if(exist('headeronly')~=1) headeronly = 0; end

fid    = fopen(fname, 'rb', 'b') ;
if(fid == -1)
  fprintf('ERROR: could not open %s for reading\n',fname);
  return;
end
v       = fread(fid, 1, 'int') ; 
if(isempty(v))
  fprintf('ERROR: problem reading fname\n');
  if(gzipped >=0) unix(sprintf('rm %s', fname)); end
end
ndim1   = fread(fid, 1, 'int') ; 
ndim2   = fread(fid, 1, 'int') ; 
ndim3   = fread(fid, 1, 'int') ; 
nframes = fread(fid, 1, 'int') ;
type    = fread(fid, 1, 'int') ; 
dof     = fread(fid, 1, 'int') ; 

if(slices(1) > 0)
  ind = find(slices > ndim3);
  if(~isempty(ind))
    fprintf('ERROR: load_mgh: some slices exceed nslices\n');
    return;
  end
end

if(frames(1) > 0)
  ind = find(frames > nframes);
  if(~isempty(ind))
    fprintf('ERROR: load_mgh: some frames exceed nframes\n');
    return;
  end
end

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
ras_good_flag = fread(fid, 1, 'short') ; 
if (ras_good_flag)
  delta  = fread(fid, 3, 'float32') ; 
  Mdc    = fread(fid, 9, 'float32') ; 
  Mdc    = reshape(Mdc,[3 3]);
  Pxyz_c = fread(fid, 3, 'float32') ; 

  D = diag(delta);

  Pcrs_c = [ndim1/2 ndim2/2 ndim3/2]'; % Should this be kept?

  Pxyz_0 = Pxyz_c - Mdc*D*Pcrs_c;

  M = [Mdc*D Pxyz_0;  ...
	0 0 0 1];
  ras_xform = [Mdc Pxyz_c; ...
	0 0 0 1];
  unused_space_size = unused_space_size - USED_SPACE_SIZE ;
end

fseek(fid, unused_space_size, 'cof') ;
nv = ndim1 * ndim2 * ndim3 * nframes;  
volsz = [ndim1 ndim2 ndim3 nframes];

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;

% Determine number of bytes per voxel
switch type
 case MRI_FLOAT,
  nbytespervox = 4;
 case MRI_UCHAR,
  nbytespervox = 1;
 case MRI_SHORT,
  nbytespervox = 2;
 case MRI_INT,
  nbytespervox = 4;
end

if(headeronly)
  fseek(fid,nv*nbytespervox,'cof');
  if(~feof(fid))
    [mr_parms count] = fread(fid,4,'float32');
    if(count ~= 4) 
      fprintf('WARNING: error reading MR params\n');
    end
  end
  fclose(fid);
  if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
  return;
end


%------------------ Read in the entire volume ----------------%
if(slices(1) <= 0 & frames(1) <= 0)
  switch type
   case MRI_FLOAT,
    vol = fread(fid, nv, 'float32') ; 
   case MRI_UCHAR,
    vol = fread(fid, nv, 'uchar') ; 
   case MRI_SHORT,
    vol = fread(fid, nv, 'short') ; 
   case MRI_INT,
    vol = fread(fid, nv, 'int') ; 
  end

  if(~feof(fid))
    [mr_parms count] = fread(fid,4,'float32');
    if(count ~= 4) 
      fprintf('WARNING: error reading MR params\n');
    end
  end
  fclose(fid) ;
  if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
  
  nread = prod(size(vol));
  if(nread ~= nv)
    fprintf('ERROR: tried to read %d, actually read %d\n',nv,nread);
    vol = [];
    return;
  end
  vol = reshape(vol,[ndim1 ndim2 ndim3 nframes]);

  return;
end

%----- only gets here if a subest of slices/frames are to be loaded ---------%


if(frames(1) <= 0) frames = [1:nframes]; end
if(slices(1) <= 0) slices = [1:ndim3]; end

nvslice = ndim1 * ndim2;
nvvol   = ndim1 * ndim2 * ndim3;
filepos0 = ftell(fid);
vol = zeros(ndim1,ndim2,length(slices),length(frames));
nthframe = 1;
for frame = frames

  nthslice = 1;
  for slice = slices
    filepos = ((frame-1)*nvvol + (slice-1)*nvslice)*nbytespervox + filepos0;
    fseek(fid,filepos,'bof');
    
    switch type
     case MRI_FLOAT,
      [tmpslice nread]  = fread(fid, nvslice, 'float32') ; 
     case MRI_UCHAR,
      [tmpslice nread]  = fread(fid, nvslice, 'uchar') ; 
     case MRI_SHORT,
      [tmpslice nread]  = fread(fid, nvslice, 'short') ; 
     case MRI_INT,
      [tmpslice nread]  = fread(fid, nvslice, 'int') ; 
    end

    if(nread ~= nvslice)
      fprintf('ERROR: load_mgh: reading slice %d, frame %d\n',slice,frame);
      fprintf('  tried to read %d, actually read %d\n',nvslice,nread);
      fclose(fid);
      if(gzipped >=0) unix(sprintf('rm %s', fname)); end
      return;
    end

    vol(:,:,nthslice,nthframe) = reshape(tmpslice,[ndim1 ndim2]);
    nthslice = nthslice + 1;
  end

  nthframe = nthframe + 1;
end

% seek to just beyond the last slice/frame %
filepos = (nframes*nvvol)*nbytespervox + filepos0;
fseek(fid,filepos,'bof');

if(~feof(fid))
  [mr_parms count] = fread(fid,5,'float32');
  if(count < 4) 
    fprintf('WARNING: error reading MR params\n');
  end
end

fclose(fid) ;
if(gzipped >=0) unix(sprintf('rm %s', fname)); end

return;
