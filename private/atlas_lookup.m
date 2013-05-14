function [label] = atlas_lookup(atlas, pos, varargin)

% ATLAS_LOOKUP determines the anatomical label of a location in the given atlas.
%
% Use as
%   atlas = atlas_init(filename);
%   label = atlas_lookup(atlas, pos, ...);
%
% Optinal input arguments should come in key-value pairs and can include
%   'queryrange'    = number, should be 1, 3, 5, 7, 9 or 11 (default = 3)
%   'inputcoord'   = 'mni' or 'tal' (default = [])
% 
% Dependent on the input coordinates and the coordinates of the atlas, the
% input positions are transformed betweem MNI and Talairach-Tournoux coordinates.
% See http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml for more details.
%
% See also ATLAS_INIT, ATLAS_MASK

% Copyright (C) 2005-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
queryrange  = ft_getopt(varargin, 'queryrange', 3);
inputcoord  = ft_getopt(varargin, 'inputcoord');

if isempty(inputcoord)
  error('you must specify inputcoord');
end

if isempty(intersect(queryrange, [1 3 5 7 9 11]))
  error('incorrect query range, should be one of [1 3 5 7 9 11]');
end

if size(pos,1)==3 && size(pos,2)~=3
  % transpose the input positions to get Nx3
  pos = pos';
end

% convert between MNI head coordinates and TAL head coordinates
% coordinates should be expressed compatible with the atlas
if     strcmp(inputcoord, 'mni') && strcmp(atlas.coord, 'tal')
  pos = mni2tal(pos')'; % this function likes 3xN 
elseif strcmp(inputcoord, 'mni') && strcmp(atlas.coord, 'mni')
  % nothing to do
elseif strcmp(inputcoord, 'tal') && strcmp(atlas.coord, 'tal')
  % nothing to do
elseif strcmp(inputcoord, 'tal') && strcmp(atlas.coord, 'mni')
  pos = tal2mni(pos')'; % this function likes 3xN 
end

num = size(pos,1);
sel = [];

% convert the atlas head coordinates into voxel coordinates
vox  = inv(atlas.transform) * [pos ones(num,1)]';
vox  = vox(1:3,:)';

for i=1:num

  % this is the center voxel
  ijk_center = vox(i,:);

  for di=(-(queryrange-1)/2):1:((queryrange-1)/2)
    for dj=(-(queryrange-1)/2):1:((queryrange-1)/2)
      for dk=(-(queryrange-1)/2):1:((queryrange-1)/2)
        
        % search in a cube around the center voxel
        ijk = round(ijk_center + [di dj dk]);

        if ijk(1)>=1 && ijk(1)<=atlas.dim(1) && ...
            ijk(2)>=1 && ijk(2)<=atlas.dim(2) && ...
            ijk(3)>=1 && ijk(3)<=atlas.dim(3)
          brick0_val = atlas.brick0(ijk(1), ijk(2), ijk(3));
          brick1_val = atlas.brick1(ijk(1), ijk(2), ijk(3));
          
          toAdd = find(atlas.descr.brick==0 & atlas.descr.value==brick0_val);
          if ~isempty(toAdd)
            sel = [sel; toAdd];
          end
          
          toAdd = find(atlas.descr.brick==1 & atlas.descr.value==brick1_val);
          if ~isempty(toAdd)
            sel = [sel; toAdd];
          end
        else
          warning('location is outside atlas volume');
        end

      end % dk
    end % dj
  end % di

end

label = unique(atlas.descr.name(sel));

