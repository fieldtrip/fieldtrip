function [vol] = yokogawa2vol(hdr)

% YOKOGAWA2VOL converts a spherical volume conductor model that can
% be present in the header of a datafile into a structure that can
% be used by FieldTrip.

% Copyright (C) 2005, Robert Oostenveld
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

% hdr = read_yokogawa_header(filename);
hdr = hdr.orig; % use the original Yokogawa header, not the FieldTrip header

if ft_hastoolbox('yokogawa_meg_reader') && (hdr.coregist.model.type == 1)
  % single sphere volume conduction model
  vol.r = hdr.coregist.model.radius;
  vol.o = [hdr.coregist.model.cx hdr.coregist.model.cy hdr.coregist.model.cz]; 
elseif isfield(hdr.mri_info, 'model_name') && strcmp(hdr.mri_info.model_name, 'SphericalModel')
  % single sphere volume conduction model
  vol.r = hdr.mri_info.r;
  vol.o = [hdr.mri_info.cx hdr.mri_info.cy hdr.mri_info.cz];
else
  error('unsupported volume conductor model');
end
