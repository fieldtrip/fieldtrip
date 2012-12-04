function lf = leadfield_interpolate(pos, vol)

% LEADFIELD_INTERPOLATE interpolates the leadfield for a source at
% an arbitrary location given the pre-computed leadfields on a regular
% grid.
%
% Use as
%   lf = leadfield_interpolate(pos, vol)

% Copyright (C) 2012, whoever works on this
%
% $Id$

% FIXME this still needs to be implemented

% compute the dipole positions in the volume conductor
xgrid = 1:vol.dim(1);
ygrid = 1:vol.dim(2);
zgrid = 1:vol.dim(3);
[x y z] = ndgrid(xgrid, ygrid, zgrid);
gridpos = warp_apply(vol.transform, [x(:) y(:) z(:)]);

gridpos(:,1) = gridpos(:,1) - pos(1);
gridpos(:,2) = gridpos(:,2) - pos(2);
gridpos(:,3) = gridpos(:,3) - pos(3);

dist = sqrt(sum(gridpos.^2, 2));
[val, indx] = min(dist);
[i1, i2, i3] = ind2sub(vol.dim, indx);

% FIXME this should of course be replaced by a better interpolation method
warning('using nearest neighbour at grid position %d, location [%f %f %f]\n', indx, gridpos(indx,:)+pos);

lf = nan(length(vol.sens.label), 3);
for i=1:length(vol.sens.label)
  % the leadfield nifti files have been mapped into memory by FT_PREPARE_VOL_SENS
  lf(i,:) = squeeze(vol.chan{i}.dat(i1, i2, i3, :));
end

