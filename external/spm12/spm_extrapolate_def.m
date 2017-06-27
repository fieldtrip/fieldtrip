function Y = spm_extrapolate_def(Y,M)
% Fill in non-finite values in a deformation field
% FORMAT Y = spm_extrapolate_def(Y,M)
% Y - the deformation field
% M - voxel-to-world transform associated with the deformation
%     (for deriving voxel sizes)
%
% This function is typically used after generating an inverse deformation,
% as these may contain missing locations.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_extrapolate_def.m 6137 2014-08-19 12:43:11Z john $

msk = isfinite(Y(:,:,:,1)); % Identify existing points
if ~all(msk(:)), % See if needed
    if nargin>=2,
        vx     = sqrt(sum(M(1:3,1:3).^2)); % Voxel sizes
    else
        vx     = [1 1 1];
    end

    % Determine affine transform to factor out
    [x1,x2,x3] = ndgrid(single(1:size(Y,1)),single(1:size(Y,2)),single(1:size(Y,3)));
    X          = cat(4,x1,x2,x3);
    M1         = spm_get_closest_affine(X,Y);
    clear X

    bnd        = spm_field('boundary'); % For tidying up afterwards
    spm_field('boundary',1);            % Free boundary conditions
    for d=1:3,
        x          = M1(d,1)*x1 + M1(d,2)*x2 + M1(d,3)*x3 + M1(d,4);
        u          = Y(:,:,:,d) - x; % Displacement field
        u(~msk)    = 0;
        u          = spm_field(single(msk),u,[vx 0 0.00001 0  1 1]); % Extrapolate displacements
        y          = Y(:,:,:,d);
        y(~msk)    = u(~msk) + x(~msk);
        Y(:,:,:,d) = y;
        clear x y u
    end
    spm_field('boundary',bnd);
end

