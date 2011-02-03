function r = save_mgh(vol, fname, M, mr_parms);
%
% save_mgh(vol,fname, M, <mr_parms>);
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz = M*[i1 i2 i3 1] where the
% indicies are 0-based
%
% mr_parms = [tr flipangle te ti]
%
% See also: load_mgh, vox2ras_0to1
%
%


%
% save_mgh.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision$
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

r = 1;

if(nargin < 2 | nargin > 4)
  msg = 'USAGE: save_mgh2(vol,fname,M)';
  return;
end

if(exist('mr_parms')~=1) mr_parms = []; end
if(isempty(mr_parms))   mr_parms = [0 0 0 0]; end
if(length(mr_parms) < 4)
  fprintf('ERROR: mr_parms length = %d, must be 4 or 5\n', ...
	  length(mr_parms));
  return;
end

% These dont appear to be used %
MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;
MRI_TENSOR = 6 ;

fid = fopen(fname, 'wb', 'b') ;
if(fid == -1)
  fprintf('ERROR: could not open %s for writing\n',fname);
  return;
end


[ndim1,ndim2,ndim3,frames] = size(vol) ;
fwrite(fid, 1, 'int') ;		% magic #
fwrite(fid, ndim1, 'int') ; 
fwrite(fid, ndim2, 'int') ; 
fwrite(fid, ndim3, 'int') ; 
fwrite(fid, frames, 'int') ;	% # of frames
if(ndims(vol) == 5)
  is_tensor = 1 ;
  fwrite(fid, MRI_TENSOR, 'int') ; % type = MRI_TENSOR
else
  is_tensor = 0 ;
  fwrite(fid, MRI_FLOAT, 'int') ;  % type = MRI_FLOAT
end

%%?????????????%%%
fwrite(fid, 1, 'int') ;          % dof (not used)
dof = fread(fid, 1, 'int') ; 

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

MdcD = M(1:3,1:3);
delta = sqrt(sum(MdcD.^2));

Mdc = MdcD./repmat(delta,[3 1]);
Pcrs_c = [ndim1/2 ndim2/2 ndim3/2 1]'; %'
Pxyz_c = M*Pcrs_c;
Pxyz_c = Pxyz_c(1:3);

fwrite(fid, 1,      'short') ;       % ras_good_flag = 1
fwrite(fid, delta,  'float32') ; 
fwrite(fid, Mdc,    'float32') ; 
fwrite(fid, Pxyz_c, 'float32') ; 

unused_space_size = UNUSED_SPACE_SIZE-2 ;
unused_space_size = unused_space_size - USED_SPACE_SIZE ;
fwrite(fid, zeros(unused_space_size,1), 'char') ;

fwrite(fid,vol,'float32');

fwrite(fid, mr_parms, 'float32') ; 
fclose(fid) ;

r = 0;

if (strcmpi(fname((length(fname)-3):length(fname)), '.MGZ') | ...
		strcmpi(fname((length(fname)-3):length(fname)), '.GZ'))

	gzipped =  round(rand(1)*10000000);
	ind = findstr(fname, '.');
	new_fname = sprintf('/tmp/tmp%d.mgh', gzipped);
	unix(sprintf('mv %s %s ; gzip %s ; mv %s.gz %s', fname, new_fname, new_fname, new_fname, fname)) ;
	fname = new_fname ;
end	
return;

