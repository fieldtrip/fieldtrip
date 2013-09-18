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

pos = warp_apply(inv(vol.transform), pos);

lf = nan(length(vol.sens.label), 3, size(pos, 1));

use_splines = (vol.chan{1}.dat.dim(end) == 6);

for i=1:length(vol.sens.label)
    for j = 1:3
        if use_splines
            lf(i, j, :) = spm_bsplins(squeeze(vol.chan{i}.dat(:, :, :, 3+j)), ...
                pos(:, 1), pos(:, 2), pos(:, 3), [4 4 4 0 0 0]);
        else
            lf(i, j, :) = spm_sample_vol(spm_vol([vol.filename{i} ',' num2str(j)]), ...
                pos(:, 1), pos(:, 2), pos(:, 3), 4);
        end
    end
end

if size(lf, 3) > 1
    lf = num2cell(lf, [size(lf, 1), 1]);
end