function hdr = load_analyze(imgfile,hdronly)
% hdr = load_analyze(imgfile,hdronly)
%
% Loads analyze volume and header. The volume is stored
% in hdr.vol. It is not permuted to swap columns and rows.
%
% hdr.dime.dim(1) = number dims
% hdr.dime.dim(2) = number of items in first dim
% hdr.dime.dim(3) = number of items in second dim
% ...
%
% hdr.dime.datatype:
%  0    unk (???) -- not supported
%  1    1 bit/pix -- not supported
%  2    8 bits
%  4   16 bits
%  8   32 bits (signed int)
% 16   32 bits (floating pt)
% 32   64 bits (2 floats) ??  -- not supported
% 64   64 bits (double)
%
% hdr.dime.pixdim(1) = physical size of first dim (eg, 3.125 mm or 2000 ms)
% hdr.dime.pixdim(2) = ...
% 
% hdr.vox2ras is the vox2ras matrix if the hdrfile is accompanied
%   by a .mat file.
%
% See also: load_analyze_hdr.m
%


%
% load_analyze.m
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

if(nargin < 1 | nargin >2)
  fprintf('hdr = load_analyze(imgfile,<hdronly>)\n');
  return;
end

if(~exist('hdronly','var')) hdronly = []; end
if(isempty(hdronly)) hdronly = 0; end

% Determine how the file name was passed:
flen = length(imgfile);
if(flen < 5)
  % Cannot have an extension, must be a base
  basename = imgfile;
else
  % Might be a file name with extension
  ext = imgfile(end-3:end);
  if(strcmp(ext,'.img') | strcmp(ext,'.hdr'))
    basename = imgfile(1:end-4);
  else
    basename = imgfile;
  end
end

imgfile0 = imgfile;
imgfile = sprintf('%s.img',basename);
hdrfile = sprintf('%s.hdr',basename);

hdr = load_analyze_hdr(hdrfile);
if(isempty(hdr)) return; end

% If only header is desired, return now
if(hdronly) return; end


% Opening img:  big or little does not seem to matter
fp = fopen(imgfile,'r',hdr.endian);
if(fp == -1) 
  fprintf('ERROR: could not open %s\n',imgfile);
  hdr = [];
  return;
end


switch(hdr.dime.datatype)
 % Note: 'char' seems to work upto matlab 7.1, but 'uchar' needed
 % for 7.2 and higher. 
 case  2, [hdr.vol nitemsread] = fread(fp,inf,'uchar');
 case  4, [hdr.vol nitemsread] = fread(fp,inf,'short');
 case  8, [hdr.vol nitemsread] = fread(fp,inf,'int');
 case 16, [hdr.vol nitemsread] = fread(fp,inf,'float');
 case 64, [hdr.vol nitemsread] = fread(fp,inf,'double');
 otherwise,
  fprintf('ERROR: data type %d not supported',hdr.dime.datatype);
  hdr = [];
  return;
end

fclose(fp);

% Get total number of voxels
dim = hdr.dime.dim(2:end);
ind0 = find(dim==0);
dim(ind0) = 1;
nvoxels = prod(dim);

% Check that that many voxels were read in
if(nitemsread ~= nvoxels) 
  fprintf('ERROR: %s, read in %d voxels, expected %d\n',...
	  imgfile,nitemsread,nvoxels);
  hdr = [];
  return;
end

hdr.vol = reshape(hdr.vol, dim');

if(hdr.dime.roi_scale ~= 0 & hdr.dime.roi_scale ~= 1)
  fprintf('fast_ldanalyze: scaling by %g\n',hdr.dime.roi_scale);
  hdr.vol = hdr.vol * hdr.dime.roi_scale;
end



return;





