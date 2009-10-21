function [vol] = read_ctf_hdm(filename);

% READ_CTF_HDM reads the head volume conductor model from a *.hdm file
%
% vol = read_ctf_hdm(filename)

% Copyright (C) 2003, Robert Oostenveld
% 
% $Log: read_ctf_hdm.m,v $
% Revision 1.2  2009/03/26 15:00:38  roboos
% remember all original details, required for 3rd gradient multisphere models
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.3  2006/02/09 08:37:23  roboos
% remove the fields SEARCH_RADIUS and HEADSHAPE_FILE regardless of where they are, do not assume that they are at location 1 and 2 (thanks to Tom)
%
% Revision 1.2  2004/06/28 07:32:55  roberto
% added units=cm for the volume model
%
% Revision 1.1  2003/03/24 12:30:42  roberto
% new implementation
%

vol   = [];
ascii = read_ctf_ascii(filename);

% remember all original details
vol.orig = ascii;

if isfield(ascii, 'MultiSphere_Data')
  chans = fieldnames(ascii.MultiSphere_Data);
  % remove the fields SEARCH_RADIUS and HEADSHAPE_FILE
  chans = chans(~(strcmp(chans, 'SEARCH_RADIUS') | strcmp(chans, 'HEADSHAPE_FILE')));
  for i=1:length(chans)
    tmp = getfield(ascii.MultiSphere_Data, chans{i});
    vol.label{i} = chans{i};
    vol.r(i) = tmp(4);
    vol.o(i, :) = tmp(1:3);
  end
  vol.r = vol.r(:);	% ensure column vector
elseif isfield(ascii, 'MEG_Sphere')
  vol.r    = ascii.MEG_Sphere.RADIUS;
  vol.o(1) = ascii.MEG_Sphere.ORIGIN_X;
  vol.o(2) = ascii.MEG_Sphere.ORIGIN_Y;
  vol.o(3) = ascii.MEG_Sphere.ORIGIN_Z;
else
  error('no headmodel information found');
end

% add the fiducials, these are in raw MRI coordinates
if isfield(ascii, 'Fid_Points')
  vol.mri.nas = ascii.Fid_Points.NASION;
  vol.mri.lpa = ascii.Fid_Points.LEFT_EAR;
  vol.mri.rpa = ascii.Fid_Points.RIGHT_EAR;
end

% add the units in which the volume conductor is defined
vol.unit = 'cm';

