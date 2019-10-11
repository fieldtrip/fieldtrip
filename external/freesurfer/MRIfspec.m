function [fspec, fstem, fmt] = MRIfspec(fstring,checkdisk)
% [fspec fstem fmt] = MRIfspec(fstring,<checkdisk>)
%
% Determine the file specification (fspec), stem (fstem), and format
% (fmt) from the given file string.
%
% A specification is the name of a file as it could exist on disk. The
% specification has the form fstem.fmt, where fmt can be mgh, mgz,
% nii, nii.gz, img, bhdr. fstring can be either an fspec or fstem. 
%
% If fstring is an fspec, then the format is determined from the
% extension.
%
% If fstring is an fstem and checkdisk=1 or NOT present, then the
% format is determined by finding a file on disk named fstem.fmt. The
% formats are searched in the following order: mgh, mgz, bhdr, img,
% nii, nii.gz. The search order is only important when there are
% multiple files that would have met the criteria, then only the first
% one is chosen. If no such file is found, then empty strings are
% returned.
%


%
% MRIfspec.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.8 $
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


fspec = [];
fstem = [];
fmt   = [];

if(nargin < 1 | nargin > 2)
  fprintf('[fspec fstem fmt] = MRIfspec(fstring,<checkdisk>)\n');
  return;
end

% If checkdisk is not passed, then do a check disk
if(~exist('checkdisk','var')) checkdisk = 1; end

inddot = max(findstr(fstring,'.'));
if(isempty(inddot))
  ext = ' ';
else
  ext = fstring(inddot+1:end);
end
    
switch(ext)
 case 'mgh',
  fspec = fstring;
  fmt = 'mgh';
  fstem = fstring(1:end-4);
  return;
 case 'mgz',
  fspec = fstring;
  fmt = 'mgz';
  fstem = fstring(1:end-4);
  return;
 case 'nii',
  fspec = fstring;
  fmt = 'nii';
  fstem = fstring(1:end-4);
  return;
 case 'gz',
  ind = findstr(fstring,'nii.gz');
  if(~isempty(ind))
    fspec = fstring;
    fmt = 'nii.gz';
    fstem = fstring(1:ind-2);
    return;
  end
 case 'img',
  fspec = fstring;
  fmt = 'img';
  fstem = fstring(1:end-4);
  return;
 case 'hdr',
  fspec = fstring;
  fmt = 'img';
  fstem = fstring(1:end-4);
  return;
 case 'bhdr',
  fspec = fstring;
  fmt = 'bhdr';
  fstem = fstring(1:end-5);
  return;
end

% If it gets here, then it cannot determine the format from an
% extension, so fstring could be a stem, so see if a file with
% stem.fmt exists. Order is imporant here.

% Do this only if checkdisk is on
if(~checkdisk) return; end

fstem = fstring;

fmt = 'mgh';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'mgz';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'bhdr';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'img';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'nii';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'nii.gz';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

% If it gets here, then could not determine format
fstem = [];
fmt = [];
fspec = [];

return;
