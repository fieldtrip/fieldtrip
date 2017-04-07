function [sts, str] = spm_check_orientations(V, verbose)
% Check the dimensions and orientations of the images
% FORMAT [sts, str] = spm_check_orientations(V [,verbose])
% V       - a struct array as returned by spm_vol
% verbose - [Default: true]
%
% sts     - status (true means OK)
% str     - string describing status, empty if OK
%
% When used without LHS, this function throws an error accordingly.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_orientations.m 5097 2012-12-06 16:08:16Z guillaume $

sts = true;
str = '';

if nargin < 2, verbose = true; end

if numel(V) <= 1 || spm_mesh_detect(V), return; end

dims = cat(1,V.dim);
if any(any(diff(dims,1,1),1))
    sts = false;
    str = 'The images do not all have the same dimensions.';
    if verbose
        fprintf('\n    ** %s **\n',strtrim(str(end,:)));
        fprintf('The function assumes that a voxel in one image corresponds with\n');
        fprintf('the same  voxel in another.   This  is not a safe assumption if\n');
        fprintf('the  image dimensions differ.   Please  ensure  that  you  have\n');
        fprintf('processed all the image data in the same way (eg. check spatial\n');
        fprintf('normalisation bounding-boxes, voxel-sizes etc).\n');
        fprintf('Here are the dimensions of the image volumes.  This list can be\n');
        fprintf('used to determine which file(s) are causing the problem.\n\n');
        for i=1:numel(V)
            fprintf('[%d %d %d]  %s\n',V(i).dim, V(i).fname);
        end
        fprintf('\n');
    end
    if ~nargout, error('The dimensions must be identical for this procedure.'); end
end

matx = reshape(cat(3,V.mat),[16,numel(V)]);
if any(any(abs(diff(matx,1,2))>1e-4))
    sts = false;
    str = char(str,'The images do not all have same orientation and/or voxel sizes.');
    if verbose
        fprintf('\n** %s **\n',strtrim(str(end,:)));
        fprintf('The function assumes that a voxel in one image  corresponds exactly\n');
        fprintf('with  the same voxel in another.   This is not a safe assumption if\n');
        fprintf('the orientation information  in the headers or .mat files says that\n');
        fprintf('the images are oriented differently. Please ensure that you process\n');
        fprintf('all data correctly. For example, you may have realigned the images,\n');
        fprintf('but not actually resliced them to be in voxel-wise alignment.\n');
        fprintf('Here are the orientation matrices of the image volumes.   This list\n');
        fprintf('can be used to determine which file(s) are causing the problem.\n\n');
        for i=1:numel(V)
            fprintf('[%g %g %g %g; %g %g %g %g; %g %g %g %g]  %s\n',...
                V(i).mat(1:3,:)', V(i).fname);
        end
        fprintf('\n');
    end
    if ~nargout, error('The orientations etc must be identical for this procedure.'); end
end
