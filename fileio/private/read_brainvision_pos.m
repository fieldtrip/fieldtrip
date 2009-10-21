function [elec] = read_brainvision_pos(filename);

% READ_BRAINVISION_POS reads electrode positions measured with the Polhemus
% tracker in one of the F.C. Donders EEG labs. The polhemus software is actually 
% not from Brainvision.
%
% Use as:
%   [elec] = read_brainvision_pos(filename)
%
% This returns an electrode structure with
%   elec.label     cell-array with electrode labels (strings)
%   elec.pnt       position of each electrode

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: read_brainvision_pos.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.1  2004/03/18 09:34:34  roberto
% also read in the fiducials if present at the end of the file
%

fid = fopen(filename, 'rt');
line = fgetl(fid);
Nchan = str2double(line);

for i=1:Nchan
  line = fgetl(fid);
  [t, r] = strtok(line);
  elec.label{i} = char(t);
  elec.pnt(i,:) = sscanf(r, '%f')';
end
elec.label = elec.label(:);

try
  % read the fiducials
  line = fgetl(fid);
  [t, r] = strtok(line);
  fiducial.label{1} = char(t);
  fiducial.pnt(1,:) = sscanf(r, '%f')';
  line = fgetl(fid);
  [t, r] = strtok(line);
  fiducial.label{2} = char(t);
  fiducial.pnt(2,:) = sscanf(r, '%f')';
  line = fgetl(fid);
  [t, r] = strtok(line);
  fiducial.label{3} = char(t);
  fiducial.pnt(3,:) = sscanf(r, '%f')';
  % add the fiducials to the electrode array
  elec.label = cat(1, elec.label(:), fiducial.label(:));
  elec.pnt   = cat(1, elec.pnt, fiducial.pnt);
catch
  % do nothing
end

fclose(fid);