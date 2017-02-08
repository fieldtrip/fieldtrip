function data = fixcoordsys(data)

% FIXCOORDSYS ensures that the coordinate system is consistently
% described. E.g. SPM and MNI are technically the same coordinate
% system, but the strings 'spm' and 'mni' are different.
%
% See also FT_DETERMINE_COORDSYS

% Copyright (C) 2017, Robert Oostenveld
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

data.coordsys = lower(data.coordsys);

% see also http://www.fieldtriptoolbox.org/faq/how_are_the_different_head_and_mri_coordinate_systems_defined
if strcmpi(data.coordsys, 'mni') || strcmpi(data.coordsys, 'spm')
  data.coordsys = 'mni';
elseif strcmpi(data.coordsys, 'ctf') || strcmpi(data.coordsys, '4d') || strcmpi(data.coordsys, 'bti')
  data.coordsys = 'ctf';
end

