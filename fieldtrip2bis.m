function fieldtrip2bis(filename, elec, mrifile)

% FIELDTRIP2BIS writes BioImage Suite .mgrid files with eletrode 
% positions in 'xyz' coordinates using a elec datatype structure and the 
% corresponding MRI volume 
%
% Use as
%   fieldtrip2bis('Subject_grid.mgrid', elec, 'Subject_MR.nii')
%
% See also BIS2FIELDTRIP, FT_WRITE_SENS, WRITE_BIOIMAGE_MGRID

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

% convert mri head coordinates to ijk coordinates
mri = ft_read_mri(mrifile);
elec_ijk.elecpos = ft_warp_apply(inv(mri.transform), elec.elecpos);

% adjust for bioimage suite indexing first voxel at [0 0 0] instead of [1 1 1]
elec_ijk.elecpos = elec_ijk.elecpos-1;

% convert ijk coordinates to xyz coordinates
elec_xyz = keepfields(mri, {'unit', 'coordsys'});
elec_xyz.label = elec.label;
elec_xyz.elecpos(:, 1) = elec_ijk.elecpos(:, 1)*mri.hdr.xsize;
elec_xyz.elecpos(:, 2) = elec_ijk.elecpos(:, 2)*mri.hdr.ysize;
elec_xyz.elecpos(:, 3) = elec_ijk.elecpos(:, 3)*mri.hdr.zsize;

% export to .mgrid file
ft_write_sens(filename, elec_xyz, 'format', 'bioimage_mgrid');
