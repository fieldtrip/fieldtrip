function elec = bis2fieldtrip(mgridfile, mrifile)

% BIS2FIELDTRIP reads BioImage Suite .mgrid files and converts them 
% into a FieldTrip-compatible elec datatype structure and converts electrode
% positions from BioImage Suite mgrid that are in 'xyz' to head coordinates
% of the corresponding MRI volume 
%
% Use as
%   elec = bis2fieldtrip('Subject_grid.mgrid', 'Subject_MR.nii')
%
% See also FIELDTRIP2BIS, FT_READ_SENS, READ_BIOIMAGE_MGRID

% Copyright (C) 2017, Arjen Stolk & Sandon Griffin
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

% import .mgrid electrode positions
elec_xyz = ft_read_sens(mgridfile, 'senstype', 'eeg');

% convert xyz coordinates to ijk coordinates
mri = ft_read_mri(mrifile);
elec_ijk.elecpos(:, 1) = elec_xyz.elecpos(:, 1)/mri.hdr.xsize;
elec_ijk.elecpos(:, 2) = elec_xyz.elecpos(:, 2)/mri.hdr.ysize;
elec_ijk.elecpos(:, 3) = elec_xyz.elecpos(:, 3)/mri.hdr.zsize;

% adjust for bioimage suite indexing first voxel at [0 0 0] instead of [1 1 1]
elec_ijk.elecpos = elec_ijk.elecpos+1;

% convert ijk coordinates to mri head coordinates
elec = keepfields(mri, {'unit', 'coordsys'});
elec.label   = elec_xyz.label;
elec.elecpos = ft_warp_apply(mri.transform, elec_ijk.elecpos);

% ensure that the elec description is up-to-date
elec = ft_datatype_sens(elec);
