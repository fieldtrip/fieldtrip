function hdr = load_nifti(niftifile,hdronly)
% hdr = load_nifti(niftifile,hdronly)
%
% Loads nifti header and volume. The volume is stored
% in hdr.vol. Columns and rows are not swapped.
%
% Handles compressed nifti (nii.gz) by issuing a unix command to
% uncompress the file to a temporary file, which is then deleted.
%
% Dimensions are in mm and msec
% hdr.pixdim(1) = physical size of first dim (eg, 3.125 mm or 2000 ms)
% hdr.pixdim(2) = ...
% 
% The sform and qform matrices are stored in hdr.sform and hdr.qform.
%
% hdr.vox2ras is the vox2ras matrix based on sform (if valid), then
% qform.
%
% Handles data structures with more than 32k cols by looking for
% hdr.dim(2) = -1 in which case ncols = hdr.glmin. This is FreeSurfer
% specific, for handling surfaces. When the total number of spatial
% voxels equals 163842, then the volume is reshaped to
% 163842x1x1xnframes. This is for handling the 7th order icosahedron
% used by FS group analysis.
%
% See also: load_nifti_hdr.m
%


%
% load_nifti.m
%
% Original Author: Doug Greve
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

hdr = [];

if(nargin < 1 || nargin > 2)
  fprintf('hdr = load_nifti(niftifile,<hdronly>)\n');
  return;
end

if(~exist('hdronly','var')) hdronly = []; end
if(isempty(hdronly)) hdronly = 0; end

% unzip if it is compressed 
ext = niftifile((strlen(niftifile)-2):strlen(niftifile));
if(strcmpi(ext,'.gz'))
  % Need to create unique file name (harder than it looks)
  %r0 = rand('state');  rand('state', sum(100*clock));
  %gzipped =  round(rand(1)*10000000 + sum(int16(niftifile))) + round(cputime);
  %rand('state',r0);
  new_niftifile = sprintf('%s.load_nifti.m.nii', tempname(fsgettmppath));
  %fprintf('Uncompressing %s to %s\n',niftifile,new_niftifile);
  gzipped = 1;
  if(strcmp(computer,'MAC') || strcmp(computer,'MACI') || ismac)
    cmd = sprintf('gunzip -c %s > %s', niftifile, new_niftifile);
  else
    cmd = sprintf('zcat %s > %s', niftifile, new_niftifile);
  end
  [status, result] = unix(cmd);
  if(status)
    fprintf('cd %s\n',pwd);
    fprintf('%s\n',cmd);
    fprintf('ERROR: %s\n',result);
    return;
  end
  niftifile = new_niftifile ;
else
  gzipped = -1 ;
end

hdr = load_nifti_hdr(niftifile);
if(isempty(hdr)) 
  if(gzipped >=0) 
    cmd = sprintf('rm -f %s', niftifile);
    [status, result] = unix(cmd); 
    if(status)
      fprintf('cd %s\n',pwd);
      fprintf('%s\n',cmd);
      fprintf('ERROR: %s\n',result);
      return;
    end
  end
  return; 
end

% Check for ico7
nspatial = prod(hdr.dim(2:4));
IsIco7 = 0;
if(nspatial == 163842) IsIco7 = 1; end

% If only header is desired, return now
if(hdronly) 
  if(gzipped >=0) unix(sprintf('rm -f %s', niftifile)); end
  if(IsIco7)
    % Reshape
    hdr.dim(2) = 163842;
    hdr.dim(3) = 1;
    hdr.dim(4) = 1;
  end
  return; 
end

% Get total number of voxels
dim = hdr.dim(2:end);
ind0 = find(dim==0);
dim(ind0) = 1;
nvoxels = prod(dim);

% Open to read the pixel data
fp = fopen(niftifile,'r',hdr.endian);

% Get past the header
fseek(fp,round(hdr.vox_offset),'bof');

switch(hdr.datatype)
 % Note: 'char' seems to work upto matlab 7.1, but 'uchar' needed
 % for 7.2 and higher. 
 case   2, dtype = 'uchar' ;
 case   4, dtype = 'short' ;
 case   8, dtype = 'int' ;
 case  16, dtype = 'float' ;
 case  64, dtype = 'double' ;
 case 256, dtype = 'int8' ;
 case 512, dtype = 'ushort' ;
 case 768, dtype = 'uint' ;
 otherwise
   fprintf('ERROR: data type %d not supported',hdr.datatype);
   hdr = [];
   fclose(fp);
   if(gzipped >=0) 
     fprintf('Deleting temporary uncompressed file %s\n',niftifile);
     unix(sprintf('rm -f %s', niftifile)); 
   end
   return;
end

% preserve volume datatype if env var is set to 1
if(getenv('FS_PRESERVE_MATLAB_VOLTYPE') == '1')
  dtype = strcat('*', dtype) ;
end 

[hdr.vol, nitemsread] = fread(fp,inf,dtype);

fclose(fp);
if(gzipped >=0) 
  %fprintf('Deleting temporary uncompressed file %s\n',niftifile);
  unix(sprintf('rm -f %s', niftifile)); 
end

% Check that that many voxels were read in
if(nitemsread ~= nvoxels) 
  fprintf('ERROR: %s, read in %d voxels, expected %d\n',...
	  niftifile,nitemsread,nvoxels);
  hdr = [];
  return;
end

if(IsIco7)
  %fprintf('load_nifti: ico7 reshaping\n');
  hdr.dim(2) = 163842;
  hdr.dim(3) = 1;
  hdr.dim(4) = 1;
  dim = hdr.dim(2:end);  
end

hdr.vol = reshape(hdr.vol, dim');
if hdr.scl_slope~=0 && ~(hdr.scl_inter==0 && hdr.scl_slope==1)
  % Rescaling is not needed if the slope==1 and intersect==0, skipping
  % this preserves the numeric class of the data
  %fprintf('Rescaling NIFTI: slope = %g, intercept = %g\n',...
  %	  hdr.scl_slope,hdr.scl_inter);
  hdr.vol = double(hdr.vol) * hdr.scl_slope  + hdr.scl_inter;
end

return;





