function dist = maskdist(vol)
%
% dist=maskdist(vol)
%
% return the distance in each voxel towards the nearest label boundaries
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    vol: a 2D or 3D array
%
% output:
%    dist: an integer array, storing the distance, in voxel unit, towards
%          the nearest boundary between two distinct non-zero voxels, the
%          zero voxels in the domain and space outside of the array
%          are also treated as a unique non-zero label. If the goal is to
%          get the minimum distance measured from the center of the voxel,
%          one should use (dist-0.5).
%
% example:
%
%    a=ones(60,60,60);
%    a(:,:,1:10)=2;
%    a(:,:,11:20)=3;
%    im=maskdist(a);
%    imagesc(squeeze(im(:,30,:)))
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if (isempty(vol))
    error('input vol can not be empty');
end

vals = unique(vol(:));
if (length(vals) > 256)
    error('it appears that your input is a gray-scale image, you must convert it to binary or labels first');
end

newvol = ones(size(vol) + 2) * max(vals) + 1;
if (ndims(vol) == 2)
    newvol(2:end - 1, 2:end - 1) = vol;
elseif (ndims(vol) == 3)
    newvol(2:end - 1, 2:end - 1, 2:end - 1) = vol;
end

vals(end + 1) = newvol(1, 1, 1);
vals(vals == 0) = [];
newvol(newvol == 0) = newvol(1, 1, 1);

dist = ones(size(newvol)) * inf;

for i = 1:length(vals(:))
    vv = (newvol == vals(i));
    vdist = bwdist(vv);
    vdist(vdist == 0) = inf;
    dist = min(dist, vdist);
end

if (ndims(vol) == 2)
    dist = dist(2:end - 1, 2:end - 1);
elseif (ndims(vol) == 3)
    dist = dist(2:end - 1, 2:end - 1, 2:end - 1);
end
