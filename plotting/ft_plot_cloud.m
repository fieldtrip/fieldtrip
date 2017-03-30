function ft_plot_cloud(pos, val, varargin)

% FT_PLOT_CLOUD visualizes spatially sparse scalar data as spherical clouds. This is
% for example useful for spectral power on depth (sEEG) electrodes.
%
% Use as
%   ft_plot_cloud(pos, val, ...)
% where the first argument are the sensor positions and the second argument are the
% values.
%
% Optional input arguments should come in key-value pairs and can include
%   'radius'          = scalar, maximum radius of cloud (default = 4)
%   'rmin'            = scalar >= 1, minimum radius of cloud if scalerad = 'yes' (default = 1)
%   'scalerad'        = scale radius with val, can be 'yes' or 'no' (default = 'yes')
%   'ptsize'          = scalar, size of points in cloud (default = 1)
%   'ptdensity'       = scalar, density of points in cloud (default = 20)
%   'ptgradient'      = scalar, degree to which density of points in cloud changes
%                       from its center, default = .5 (uniform density)
%   'colormap'        = colormap for functional data, see COLORMAP
%   'colorgrad'       = 'white' or a scalar (e.g. 1), degree to which color of points
%                       in cloud changes from its center
%   'clim'            = 1x2 vector specifying the min and max for the colorscale
%   'unit'            = string, convert the sensor array to the specified geometrical units (default = [])
%
% See also FT_ELECTRODEPLACEMENT, FT_PLOT_TOPO, FT_PLOT_TOPO3D

% Copyright (C) 2017, Arjen Stolk, Sandon Griffin
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

% run some checks
if size(pos,2)~=3
  error('pos has to be an Nx3 array')
end

if isempty(val)
  val = ones(size(pos,1),1); % vector of ones
end

assert(isrow(val) || iscolumn(val), 'values should be represented as a single vector')
val = val(:); % ensure it is a column

% the unit of geometry is needed for the other defaults
% estimate the unit of geometry of the positions, FIXME this fails for a single or few points
posunit = ft_estimate_units(range(pos));
% determine the desired unit of geometry
unit = ft_getopt(varargin, 'unit', posunit);
% convert the sensor positions into the desired units
pos = pos * ft_scalingfactor(posunit, unit);

% get the input arguments
radius          = ft_getopt(varargin, 'radius', 4 * ft_scalingfactor('mm', unit));
rmin            = ft_getopt(varargin, 'rmin',   1 * ft_scalingfactor('mm', unit));
scalerad        = ft_getopt(varargin, 'scalerad', 'yes');
ptsize          = ft_getopt(varargin, 'ptsize', 1);
ptdens          = ft_getopt(varargin, 'ptdensity', 20);
ptgrad          = ft_getopt(varargin, 'ptgradient', .5);
cmap            = ft_getopt(varargin, 'colormap', 'parula');
cgrad           = ft_getopt(varargin, 'colorgrad', 'white');
clim            = ft_getopt(varargin, 'clim');

if rmin < 1 * ft_scalingfactor('mm', unit)
  error('cfg.rmin must be equal or larger than 1 mm');
end
if isempty(clim)
  clim = [min(val) max(val)]; % use the data
end

% functional data scaling factors
if ischar(cmap)
  if strcmp(cmap, 'default')
    cmapsc = get(0, 'DefaultFigureColormap');
  else
    cmapsc = feval(cmap, 201); % an odd number
  end
else
  cmapsc = cmap;
end

cmid    = size(cmapsc,1)/2;               % colorbar middle
colscf  = val / max(abs(clim));           % color: between -1 and 1 (used when cgrad = 'white')
colscf(colscf>1)=1; colscf(colscf<-1)=-1; % clamp values outside the [-1 1] range
radscf = abs( val / max(abs(clim)) );     % radius: between 0 and 1 (used when cgrad = a scalar)
radscf(radscf>1)=1; radscf(radscf<0)=0;   % clamp values outside the [0 1] range

% generate point cloud(s)
hold on;
for n = 1:size(pos,1) % cloud loop
  
  % point cloud with radius scaling
  rng(0, 'twister'); % random number generator
  if strcmp(scalerad, 'yes')
    rmax = rmin+(radius-rmin)*radscf(n); % maximum radius of this cloud
  else
    rmax = radius; % each cloud has the same radius
  end
  npoints   = round(ptdens*(4/3)*pi*rmax^3);  % number of points based on cloud volume
  elevation = asin(2*rand(npoints,1)-1);      % elevation angle for each point
  azimuth   = 2*pi*rand(npoints,1);           % azimuth angle for each point
  radii     = rmax*(rand(npoints,1).^ptgrad); % radius value for each point
  radii     = sort(radii);                    % sort radii in ascending order so they are plotted from inside out
  [x,y,z]   = sph2cart(azimuth, elevation, radii); % convert to Carthesian
  
  % color axis with radius scaling
  if strcmp(cgrad, 'white')                   % color runs up to white
    indx  = ceil(cmid) + sign(colscf(n))*floor(abs(colscf(n)*cmid));
    indx  = max(min(indx,size(cmapsc,1)),1);  % index should fall within the colormap
    fcol  = cmapsc(indx,:);                   % color [Nx3]
    ptcol = [linspace(fcol(1), 1, npoints)' linspace(fcol(2), 1, npoints)' linspace(fcol(3), 1, npoints)'];
  elseif isscalar(cgrad)                      % color runs down towards colorbar middle
    rnorm = radii/rmax;                       % normalized radius
    if radscf(n)>=.5                          % extreme values
      ptcol = val(n) - (flip(1-rnorm).^inv(cgrad))*val(n); % scaled fun [Nx1]
    elseif radscf(n)<.5                       % values closest to zero
      ptcol = val(n) + (flip(1-rnorm).^inv(cgrad))*abs(val(n)); % scaled fun [Nx1]
    end
  else
    error('color gradient should be either ''white'' or a scalar determining color falloff')
  end
  
  % draw the points
  scatter3(x+pos(n,1), y+pos(n,2), z+pos(n,3), ptsize, ptcol, '.');
  
end % end cloud loop

% colorbar settings
colormap(cmap);
if ~isempty(clim) && clim(2)>clim(1)
  caxis(gca, clim);
end

% axis settings
axis off
axis vis3d
axis equal
