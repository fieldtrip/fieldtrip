function [hx, hy, hz] = ft_plot_ortho(dat, varargin)

% FT_PLOT_ORTHO plots 3 orthographic slices through a 3-D volume and interpolates if needed
%
% Use as
%   ft_plot_ortho(dat, ...)
% or
%   ft_plot_ortho(dat, mask, ...)
% where dat and mask are equal-sized 3-D arrays.
%
% Additional options should be specified in key-value pairs and can be
%   'style'        = string, 'subplot' or 'intersect' (default = 'subplot')
%   'orientation'  = 3x3 matrix specifying the directions orthogonal through the planes which will be plotted
%   'parents'      = (optional) 3-element vector containing the handles of the axes for the subplots (when style = 'subplot')
%   'surfhandle'   = (optional) 3-element vector containing the handles of the surfaces for each of the sublots (when style = 'subplot'). Parents and surfhandle are mutually exclusive
%   'update'       = (optional) 3-element boolean vector with the axes that should be updated (default = [true true true])
%
% The following options are supported and passed on to FT_PLOT_SLICE
%   'clim'                = [min max], lower and upper color limits
%   'transform'           = 4x4 homogeneous transformation matrix specifying the mapping from voxel space to the coordinate system in which the data are plotted
%   'location'            = 1x3 vector specifying the intersection point at which the three slices will be plotted. The coordinates should be expressed in the coordinate system of the data. 
%   'datmask'             = 3D-matrix with the same size as the matrix dat, serving as opacitymap if the second input argument to the function contains a matrix, this will be used as the mask
%   'maskstyle'           = string, 'opacity' or 'colormix', defines the rendering
%   'background'          = needed when maskstyle is 'colormix', 3D-matrix with the same size as the data matrix, serving as grayscale image that provides the background
%   'interpmethod'        = string specifying the method for the interpolation, see INTERPN (default = 'nearest')
%   'colormap'            = string, see COLORMAP
%   'unit'                = string, can be 'm', 'cm' or 'mm (default is automatic)
%   'intersectmesh'       = triangulated mesh, see FT_PREPARE_MESH
%   'intersectcolor'      = string, color specification
%   'intersectlinestyle'  = string, line specification 
%   'intersectlinewidth'  = number
%
% See also FT_PLOT_SLICE, FT_PLOT_MONTAGE, FT_SOURCEPLOT

% Copyrights (C) 2010, Jan-Mathijs Schoffelen
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

% parse first input argument(s). it is either
% (dat, varargin)
% (dat, msk, varargin)
% (dat, [], varargin)
% this is done in ft_plot_slice

sellist = 1:numel(varargin);
if ~isempty(sellist)
  if isempty(varargin{1}) || isnumeric(varargin{1})
    sellist(1) = [];
  end
end

% get the optional input arguments
% other options such as location and transform are passed along to ft_plot_slice
style     = ft_getopt(varargin(sellist), 'style', 'subplot');
ori       = ft_getopt(varargin(sellist), 'orientation', eye(3));

if strcmp(style, 'subplot')
  parents    = ft_getopt(varargin(sellist), 'parents');
  surfhandle = ft_getopt(varargin(sellist), 'surfhandle');
  patchhandle = ft_getopt(varargin(sellist), 'patchhandle');
  update     = ft_getopt(varargin(sellist), 'update', [true true true]);
  if ~isempty(surfhandle) && ~isempty(parents)
    ft_error('if specifying handles, you should either specify handles to the axes or to the surface objects, not both');
  end
end

if ~isa(dat, 'double')
  dat = cast(dat, 'double');
end

% determine the orientation key-value pair
keys = varargin(sellist(1:2:end));

sel  = find(strcmp('orientation', keys));
if isempty(sel)
  % add orientation key-value pair if it does not exist
  sel             = numel(keys)+1;
  varargin{2*sel-1} = 'orientation';
  varargin{2*sel} = [];
end

switch style
  case 'subplot'
    if isempty(parents) && isempty(surfhandle)
      Hx = subplot(2,2,1);
      Hy = subplot(2,2,2);
      Hz = subplot(2,2,4);
    elseif ~isempty(parents) && isempty(surfhandle)
      Hx = parents(1);
      Hy = parents(2);
      Hz = parents(3);
    elseif isempty(parents) && ~isempty(surfhandle)
      % determine the parents from the surface handle and use the
      % surfhandle for efficient visualization (overwriting existing data)
      if update(1), Hx = get(surfhandle(1), 'parent'); else Hx = []; end
      if update(2), Hy = get(surfhandle(2), 'parent'); else Hy = []; end
      if update(3), Hz = get(surfhandle(3), 'parent'); else Hz = []; end
    end
    
    if ~isempty(Hx)
      if ~isempty(surfhandle) && update(1)
        varargin(sellist) = ft_setopt(varargin(sellist), 'surfhandle', surfhandle(1));
      end
      if ~isempty(patchhandle) && update(1)
        varargin(sellist) = ft_setopt(varargin(sellist), 'patchhandle', patchhandle(1));
      end
      % swap the first 2 dimensions because of meshgrid vs ndgrid issues
      varargin{2*sel} = ori(2,:);
      set(gcf,'currentaxes',Hx);
      hx = ft_plot_slice(dat, varargin{:});
      set(Hx, 'view', [0 0]); %, 'xlim', [0.5 size(dat,1)-0.5], 'zlim', [0.5 size(dat,3)-0.5]);
      if isempty(parents)
        % only change axis behavior if no parents are specified
        axis off
      end
    end
    
    if ~isempty(Hy)
      if ~isempty(surfhandle) && update(2)
        varargin(sellist) = ft_setopt(varargin(sellist), 'surfhandle', surfhandle(2));
      end
      if ~isempty(patchhandle) && update(2)
        varargin(sellist) = ft_setopt(varargin(sellist), 'patchhandle', patchhandle(2));
      end
      varargin{2*sel} = ori(1,:);
      set(gcf,'currentaxes',Hy);
      hy = ft_plot_slice(dat, varargin{:});
      set(Hy, 'view', [90 0]); %, 'ylim', [0.5 size(dat,2)-0.5], 'zlim', [0.5 size(dat,3)-0.5]);
      if isempty(parents)
        % only change axis behavior if no parents are specified
        axis off
      end
    end
    
    if ~isempty(Hz)
      if ~isempty(surfhandle) && update(3)
        varargin(sellist) = ft_setopt(varargin(sellist), 'surfhandle', surfhandle(3));
      end
      if ~isempty(patchhandle) && update(3)
        varargin(sellist) = ft_setopt(varargin(sellist), 'patchhandle', patchhandle(3));
      end
      varargin{2*sel} = ori(3,:);
      set(gcf,'currentaxes',Hz);
      hz = ft_plot_slice(dat, varargin{:});
      set(Hz, 'view', [0 90]); %, 'xlim', [0.5 size(dat,1)-0.5], 'ylim', [0.5 size(dat,2)-0.5]);
      if isempty(parents)
        % only change axis behavior if no parents are specified
        axis off
      end
    end
    
  case 'intersect'
    holdflag = ishold;
    if ~holdflag
      hold on
    end
    
    varargin{2*sel} = ori(1,:);
    hx = ft_plot_slice(dat, varargin{:});
    
    varargin{2*sel} = ori(2,:);
    hy = ft_plot_slice(dat, varargin{:});
    
    varargin{2*sel} = ori(3,:);
    hz = ft_plot_slice(dat, varargin{:});
    axis equal; axis tight; axis off;axis vis3d
    view(3);
    
    if ~holdflag
      hold off
    end
    
  otherwise
    ft_error('unsupported style %s', style);
    
end % switch style

% if strcmp(interactive, 'yes')
%   flag = 1;
%   while flag
%
%   end
% end
