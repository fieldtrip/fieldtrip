function [vol] = read_ctf_hdm(filename)

% READ_CTF_HDM reads the head volume conductor model from a *.hdm file
%
% vol = read_ctf_hdm(filename)

% Copyright (C) 2003, Robert Oostenveld
% 
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

vol   = [];
ascii = read_ctf_ascii(filename);

% remember all original details
vol.orig = ascii;

if isfield(ascii, 'MultiSphere_Data')
  chans = fieldnames(ascii.MultiSphere_Data);
  % remove the fields SEARCH_RADIUS and HEADSHAPE_FILE
  chans = chans(~strcmp(chans, 'SEARCH_RADIUS'));
  chans = chans(~strcmp(chans, 'HEADSHAPE_FILE'));
  chans = chans(~strcmp(chans, 'SURFACE_TYPE'));
  chans = chans(~strcmp(chans, 'HEADPOS'));
  for i=1:length(chans)
    tmp = getfield(ascii.MultiSphere_Data, chans{i});
    vol.label{i} = chans{i};
    vol.r(i)     = tmp(4);
    vol.o(i, :)  = tmp(1:3);
  end
  vol.r = vol.r(:); % ensure column vector
elseif isfield(ascii, 'MEG_Sphere')
  vol.r    = ascii.MEG_Sphere.RADIUS;
  vol.o(1) = ascii.MEG_Sphere.ORIGIN_X;
  vol.o(2) = ascii.MEG_Sphere.ORIGIN_Y;
  vol.o(3) = ascii.MEG_Sphere.ORIGIN_Z;
else
  ft_error('no headmodel information found');
end

% add the fiducials, these are in raw MRI coordinates
if isfield(ascii, 'Fid_Points')
  vol.mri.nas = ascii.Fid_Points.NASION;
  vol.mri.lpa = ascii.Fid_Points.LEFT_EAR;
  vol.mri.rpa = ascii.Fid_Points.RIGHT_EAR;
end

% add the units in which the volume conductor is defined
vol.unit = 'cm';

