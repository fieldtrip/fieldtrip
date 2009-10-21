function [hc] = read_ctf_hc(filename);

% READ_CTF_HC reads the MEG headcoil marker positions from an ascii file
% and computes the coordinate transformation required to get from from
% dewar to head-coordinates
%
% the definition of head coordinates is according to CTF standard:
% - the origin is exactly between LPA and RPA
% - the positive x-axis goes throught NAS
% - the positive y-axis goes (approximately) through LPA
% - the positive z-axis goes up, orthogonal to the x- and y-axes
%
% hc = read_ctf_hc(filename)
%
% returns a structure with the following fields
%   hc.dewar.nas	marker positions relative to dewar 
%   hc.dewar.lpa
%   hc.dewar.rpa
%   hc.head.nas		marker positions relative to head (measured) 
%   hc.head.lpa
%   hc.head.rpa
%   hc.standard.nas	marker positions relative to head (expected)
%   hc.standard.lpa
%   hc.standard.rpa
% and
%   hc.affine		parameter for affine transformation (1x12)
%   hc.homogenous	homogenous transformation matrix (4x4, see warp3d)
%   hc.translation	translation vector (1x3)
%   hc.rotation		rotation matrix (3x3)
% 
% Gradiometer positions can be transformed into head coordinates using the 
% homogeneous transformation matrix, or using the affine parameters and
% the warp3d function from the WARPING toolbox

% Copyright (C) 2002, Robert Oostenveld
% 
% $Log: read_ctf_hc.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.4  2005/06/01 08:04:24  roboos
% removed the conversion from cm to mm, output now is in mm for consistency with other ctf functions
%
% Revision 1.3  2003/03/12 16:15:18  roberto
% added extra warning for mm units
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

global fb

hc.standard.nas = [0 0 0];
hc.standard.lpa = [0 0 0];
hc.standard.rpa = [0 0 0];
hc.dewar.nas    = [0 0 0];
hc.dewar.lpa    = [0 0 0];
hc.dewar.rpa    = [0 0 0];
hc.head.nas     = [0 0 0];
hc.head.lpa     = [0 0 0];
hc.head.rpa     = [0 0 0];

fid = fopen(filename, 'r');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

fseek(fid, 0, 'bof');
line = [];

if fb
  fprintf('reading standard coil positions with respect to the dewar\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard coil positions with respect to the dewar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~strcmp(line, 'standard nasion coil position relative to dewar (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.nas(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.nas(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.nas(3) = str2num(r(2:end));
% NOTE THAT THERE IS AN TYPING ERROR IN SOME CTF FILES WHICH I HAVE TO REPRODUCE HERE (staNdard)
while ~(strcmp(line, 'stadard left ear coil position relative to dewar (cm):') | ...
        strcmp(line, 'standard left ear coil position relative to dewar (cm):'))
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.lpa(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.lpa(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.lpa(3) = str2num(r(2:end));
while ~strcmp(line, 'standard right ear coil position relative to dewar (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.rpa(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.rpa(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.standard.rpa(3) = str2num(r(2:end));

if fb
  fprintf('reading measured coil positions with respect to the dewar\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measured coil positions with respect to the dewar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~strcmp(line, 'measured nasion coil position relative to dewar (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.nas(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.nas(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.nas(3) = str2num(r(2:end));
while ~strcmp(line, 'measured left ear coil position relative to dewar (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.lpa(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.lpa(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.lpa(3) = str2num(r(2:end));
while ~strcmp(line, 'measured right ear coil position relative to dewar (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.rpa(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.rpa(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.dewar.rpa(3) = str2num(r(2:end));

if fb
  fprintf('reading measured coil positions with respect to the head\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measured coil positions with respect to the head
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~strcmp(line, 'measured nasion coil position relative to head (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.nas(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.nas(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.nas(3) = str2num(r(2:end));
while ~strcmp(line, 'measured left ear coil position relative to head (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.lpa(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.lpa(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.lpa(3) = str2num(r(2:end));
while ~strcmp(line, 'measured right ear coil position relative to head (cm):')
  line = fgetl(fid);
  if ~ischar(line) & line==-1, error('premature end of file'), end
end
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.rpa(1) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.rpa(2) = str2num(r(2:end));
line = fgetl(fid); [t, r] = strtok(line, '='); hc.head.rpa(3) = str2num(r(2:end));

fclose(fid);

if fb
  fprintf('computing transformation matrix\n');
end

% compute the direction of the head coordinate axes in dewar coordinates
d_x = hc.dewar.nas - (hc.dewar.lpa + hc.dewar.rpa)/2;
d_z = cross(d_x, hc.dewar.lpa - hc.dewar.rpa);
d_y = cross(d_z, d_x);
d_x = d_x / norm(d_x);
d_y = d_y / norm(d_y);
d_z = d_z / norm(d_z);

% compute the translation and rotation which are neccessary to transform
% any given location from the dewar to the head coordinate system
hc.translation = -(hc.dewar.lpa + hc.dewar.rpa)/2;
hc.rotation = inv([d_x' d_y' d_z']);

% construct the homogenous coordinate transformation matrix
% from the subsequent translation and rotation
hc.homogenous = eye(4,4);
hc.homogenous(1:3,4)   = hc.rotation * hc.translation';
hc.homogenous(1:3,1:3) = hc.rotation;

% construct a vector with the 12 parameters for an affine warp
hc.affine = hc.homogenous';
hc.affine = hc.affine(1:12);

