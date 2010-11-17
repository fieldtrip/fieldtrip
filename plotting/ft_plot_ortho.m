function ft_plot_ortho(dat, varargin)

% FT_PLOT_ORTHO plots a 3 orthographic cuts through a 3-D volume
%
% Use as
%   ft_plot_ortho(dat, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'style'        'subplot' (default) or 'intersect'
%   'transform'    a 4x4 homogeneous transformation matrix specifying the mapping from
%                    voxel space to the coordinate system in which the data are plotted.
%   'location'     a 1x3 vector specifying a point on the plane which will be plotted
%                    the coordinates are expressed in the coordinate system in which the 
%                    data will be plotted. location defines the origin of the plane 
%   'orientation'  a 3x3 matrix specifying the directions orthogonal through the planes
%                    which will be plotted.
%   'datmask'      a 3D-matrix with the same size as the matrix dat, serving as opacitymap
%   'interpmethod' a string specifying the method for the interpolation, default = 'nearest' 
%                    see INTERPN
%   'colormap'    
%
%   'interplim'
%
% This function uses FT_PLOT_SLICE 

% Copyrights (C) 2010, Jan-Mathijs Schoffelen
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

style     = keyval('style',       varargin); if isempty(style),      style      = 'subplot'; end
loc       = keyval('location',    varargin); if isempty(loc),        loc        = [0 0 0];   end
ori       = keyval('orientation', varargin); if isempty(ori),        ori        = eye(3);    end   
transform = keyval('transform',   varargin); if isempty(transform),  transform  = eye(4);    end
if ~strcmp(class(dat), 'double'),
  dat       = cast(dat, 'double');
end

clear ft_plot_slice

% add orientation key-value pair if it does not exist
keys = varargin(1:2:end);
sel  = strmatch('orientation', keys);
if isempty(sel)
  sel           = numel(varargin)+1;
  varargin{sel} = 'orientation';
  varargin{sel+1} = [];
end

switch style
case 'subplot'
  
  Hx = subplot(2,2,1);
  varargin{sel+1} = ori(1,:);
  hx = ft_plot_slice(dat, varargin{:});
  view([90 0]);
  axis equal;axis tight
  
  Hy = subplot(2,2,2);
  varargin{sel+1} = ori(2,:);
  hy = ft_plot_slice(dat, varargin{:});
  view([0 0]);
  axis equal;axis tight
  
  Hz = subplot(2,2,4);
  varargin{sel+1} = ori(3,:);
  hz = ft_plot_slice(dat, varargin{:});
  view([0 90]);
  axis equal;axis tight
  
case 'intersect'
  hold on;
  
  varargin{sel+1} = ori(1,:);
  hx = ft_plot_slice(dat, varargin{:});
  
  varargin{sel+1} = ori(2,:);
  hy = ft_plot_slice(dat, varargin{:});
  
  varargin{sel+1} = ori(3,:);
  hz = ft_plot_slice(dat, varargin{:});
  axis equal; axis tight; axis off;axis vis3d
  view(3);

otherwise
end
