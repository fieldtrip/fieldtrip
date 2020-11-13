function [label] = atlas_lookup(atlas, pos, varargin)

% ATLAS_LOOKUP determines the anatomical label of a location in the given atlas.
%
% Use as
%   label = atlas_lookup(atlas, pos, ...);
%
% Optinal input arguments should come in key-value pairs and can include
%   'method'       = 'sphere' (default) searches surrounding voxels in a sphere
%                    'cube' searches surrounding voxels in a cube
%   'queryrange'   = number, should be 1, 3, 5, 7, 9 or 11 (default = 3)
%   'coordsys'     = 'mni' or 'tal' (default = [])
%
% Dependent on the coordinates if the input points and the coordinates of the atlas,
% the input positions are transformed betweem MNI and Talairach-Tournoux coordinates.
% See http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml for more details.

% Copyright (C) 2005-2020, Robert Oostenveld
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

% get the optional input arguments
method      = ft_getopt(varargin, 'method', 'sphere');
queryrange  = ft_getopt(varargin, 'queryrange', 3);
coordsys    = ft_getopt(varargin, 'coordsys');

if isempty(coordsys)
  ft_error('you must specify coordsys');
end

if isempty(intersect(queryrange, 1:2:queryrange))
  ft_error('incorrect query range, should be an odd number');
end

if size(pos,1)==3 && size(pos,2)~=3
  % transpose the input positions to get Nx3
  pos = pos';
end

% determine which field(s) to use to look up the labels,
% and whether these are boolean or indexed
fn = fieldnames(atlas);
isboolean = false(numel(fn),1);
isindexed = false(numel(fn),1);
for i=1:length(fn)
  if islogical(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
    isboolean(i) = true;
  elseif isnumeric(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
    isindexed(i) = true;
  end
end
if any(isindexed)
  % let the indexed prevail
  fn = fn(isindexed);
  isindexed = 1;
elseif any(isboolean)
  % use the boolean
  fn = fn(isboolean);
  isindexed = 0;
end

% convert between MNI head coordinates and TAL head coordinates
% coordinates should be expressed compatible with the atlas
if     strcmp(coordsys, 'mni') && strcmp(atlas.coordsys, 'tal')
  pos = mni2tal(pos')'; % this function likes 3xN
elseif strcmp(coordsys, 'mni') && strcmp(atlas.coordsys, 'mni')
  % nothing to do
elseif strcmp(coordsys, 'tal') && strcmp(atlas.coordsys, 'tal')
  % nothing to do
elseif strcmp(coordsys, 'tal') && strcmp(atlas.coordsys, 'mni')
  pos = tal2mni(pos')'; % this function likes 3xN
elseif ~strcmp(coordsys, atlas.coordsys)
  ft_error('the mismatch between the coordinate system in the atlas and the coordinate system in the data cannot be resolved');
end

num = size(pos,1);
sel = cell(1,numel(fn));
label = {};

% convert the atlas head coordinates into voxel coordinates
vox  = ft_warp_apply(inv(atlas.transform), pos);

for i=1:num
  
  % this is the center voxel
  ijk_center = vox(i,:);
  
  if isindexed
    if strcmp(method, 'sphere')
      % search in a sphere around the center voxel
      
      % first, identify the voxels (x,y,z) in a sphere around the center voxel
      [x, y, z] = sphere(1000);
      ori = [0 0 0];
      ptswithinq = [];
      for r = 0:(queryrange/2 - 0.5)
        xs = round(r*x(:));
        ys = round(r*y(:));
        zs = round(r*z(:));
        pts = unique([xs ys zs], 'rows');
        
        spherepts = [];
        for a = 1:size(pts,1)
          d2ori = sqrt(sum((pts(a, :)-ori).^2));
          if d2ori <= r
            spherepts = [spherepts; pts(a, :)];
          end
        end
        
        ptswithinr = []; % voxels located at a radius of r voxels from the center voxel
        for n = 1:size(spherepts,1)
          ptswithinr = [ptswithinr; spherepts(n, 1)+ijk_center(1), spherepts(n, 2)+ijk_center(2), spherepts(n,3)+ijk_center(3)];
        end
        
        ptswithinq = [ptswithinq; ptswithinr]; % voxels within a radius of queryrange from the center voxel
      end
      ptswithinq = unique(round(ptswithinq), 'rows');
      
      for n = 1:size(ptswithinq,1)
        ijk = ptswithinq(n, :);
        if ijk(1)>=1 && ijk(1)<=atlas.dim(1) && ...
            ijk(2)>=1 && ijk(2)<=atlas.dim(2) && ...
            ijk(3)>=1 && ijk(3)<=atlas.dim(3)
          for k=1:numel(fn)
            sel{k} = [sel{k}; atlas.(fn{k})(ijk(1), ijk(2), ijk(3))];
          end
        else
          ft_warning('location is outside atlas volume');
        end
      end
      
    elseif strcmp(method, 'cube')
      % search in a cube around the center voxel
      for di=(-(queryrange-1)/2):1:((queryrange-1)/2)
        for dj=(-(queryrange-1)/2):1:((queryrange-1)/2)
          for dk=(-(queryrange-1)/2):1:((queryrange-1)/2)
            
            % search in a cube around the center voxel
            ijk = round(ijk_center + [di dj dk]);
            
            if ijk(1)>=1 && ijk(1)<=atlas.dim(1) && ...
                ijk(2)>=1 && ijk(2)<=atlas.dim(2) && ...
                ijk(3)>=1 && ijk(3)<=atlas.dim(3)
              for k=1:numel(fn)
                sel{k} = [sel{k}; atlas.(fn{k})(ijk(1), ijk(2), ijk(3))];
              end
              %brick0_val = atlas.brick0(ijk(1), ijk(2), ijk(3));
              %brick1_val = atlas.brick1(ijk(1), ijk(2), ijk(3));
              %sel = [sel; find(atlas.descr.brick==0 & atlas.descr.value==brick0_val)];
              %sel = [sel; find(atlas.descr.brick==1 & atlas.descr.value==brick1_val)];
            else
              ft_warning('location is outside atlas volume');
            end % k
            %FIXME the three loops can probably be easily vectorized
          end % dk
        end % dj
      end % di
    end
    
    for k = 1:numel(fn)
      if ~isempty(sel{k})
        % Get rid of zeros in sel{k}
        for t = numel(sel{k}):-1:1
          if sel{k}(t) == 0
            sel{k}(t) = [];
          end
        end
        label = [label; atlas.([fn{k} 'label'])(sel{k})]; % by using setdiff and/or unique, the count for each label is lost, and ft_volumelookup cannot provide an accurate number for labels.count
      end
    end
  else
    ft_error('support for atlases that have a probabilistic segmentationstyle is not supported yet');
  end
end


%label = unique(atlas.descr.name(sel));

