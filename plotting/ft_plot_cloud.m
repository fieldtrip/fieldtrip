function ft_plot_cloud(pos, val, varargin)

% FT_PLOT_CLOUD visualizes spatially sparse scalar data as spherical clouds and
% optionally 2D slices through the spherical clouds. This is for example useful for
% spectral power on depth (sEEG) electrodes.
%
% Use as
%   ft_plot_cloud(pos, val, ...)
% where the first argument are the sensor positions and the second argument are the
% sensor values.
%
% Optional input arguments should come in key-value pairs and can include
%   'radius'             = scalar, maximum radius of cloud (default = 4)
%   'rmin'               = scalar >= 1, minimum radius of cloud if scalerad = 'yes' (default = 1)
%   'scalerad'           = scale radius with val, can be 'yes' or 'no' (default = 'yes')
%   'colormap'           = colormap for functional data, see COLORMAP
%   'colorgrad'          = 'white' or a scalar (e.g. 1), degree to which color of points
%                          in cloud changes from its center
%   'clim'               = 1x2 vector specifying the min and max for the colorscale
%   'unit'               = string, convert the sensor array to the specified geometrical units (default = [])
%   'mesh'               = string or Nx1 cell array, triangulated mesh(es), see FT_PREPARE_MESH
%   'slice'              = requires 'mesh' as input (default = 'none')
%                          '2d', plots 2D slices through the cloud with an outline of the mesh
%                          '3d', draws an outline around the mesh at a particular slice
%   'slicetype'          = 'surf' plots the slices as a surface
%                          'point' (default) plots the slices as points
%
% The following inputs apply when 'slice' = 'none' or '3d', or '2dslicetype' = 'point'
%   'ptsize'             = scalar, size of points in cloud (default = 1)
%   'ptdensity'          = scalar, density of points in cloud (default = 20)
%   'ptgradient'         = scalar, degree to which density of points in cloud changes
%                          from its center, default = .5 (uniform density)
%
% The following inputs apply when 'slice' = '2d' or '3d'
%   'ori'                = 'x', 'y', or 'z', specifies the orthogonal plane which will be plotted (default = 'y')
%   'slicepos'           = 'auto' or Nx1 vector specifying the position of the
%                          slice plane along the orientation axis (default = 'auto': chooses slice(s) with
%                          the most data)
%   'nslices'            = scalar, number of slices to plot if 'slicepos' = 'auto (default = 1)
%   'minspace'           = scalar, minimum spacing between slices if nslices>1
%                          (default = 1)
%   'intersectcolor'     = string, Nx1 cell array, or Nx3 vector specifying line color (default = 'k')
%   'intersectlinestyle' = string or Nx1 cell array, line style specification (default = '-')
%   'intersectlinewidth' = scalar or Nx1 vector, line width specification (default = 2)
%
% See also FT_ELECTRODEPLACEMENT, FT_PLOT_TOPO, FT_PLOT_TOPO3D

% The following inputs apply when 'slicetype' = 'surf'
%   'ncirc'           = scalar, number of concentric circles to plot for each
%                       cloud slice (default = 15) make this hidden or scale
%   'scalealpha'      = 'yes' or 'no', scale the maximum alpha value of the center circle
%                       with distance from center of cloud

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

% estimate the unit of geometry of the positions (needed for the other defaults)
posunit = ft_estimate_units(range(pos).*2); % FIXME times 2 otherwise this fails for single/few points
% determine the desired unit of geometry
unit = ft_getopt(varargin, 'unit', posunit);
% convert the sensor positions into the desired units
pos = pos * ft_scalingfactor(posunit, unit);

% get the generic input arguments
radius             = ft_getopt(varargin, 'radius', 4 * ft_scalingfactor('mm', unit));
rmin               = ft_getopt(varargin, 'rmin',   1 * ft_scalingfactor('mm', unit));
scalerad           = ft_getopt(varargin, 'scalerad', 'yes');
cgrad              = ft_getopt(varargin, 'colorgrad', 'white');
if ft_platform_supports('parula')
  cmap             = ft_getopt(varargin, 'colormap', 'parula');
else
  cmap             = ft_getopt(varargin, 'colormap', 'jet');
end
clim               = ft_getopt(varargin, 'clim');
meshplot           = ft_getopt(varargin, 'mesh');

% point related inputs
ptsize             = ft_getopt(varargin, 'ptsize', 1);
ptdens             = ft_getopt(varargin, 'ptdensity', 20);
ptgrad             = ft_getopt(varargin, 'ptgradient', .5);

% slice related inputs
sli                = ft_getopt(varargin, 'slice', 'none');
slicetype          = ft_getopt(varargin, 'slicetype', 'point');
ori                = ft_getopt(varargin, 'ori', 'y');
slicepos           = ft_getopt(varargin, 'slicepos', 'auto');
nslices            = ft_getopt(varargin, 'nslices', 1);
minspace           = ft_getopt(varargin, 'minspace', 1);
intersectcolor     = ft_getopt(varargin, 'intersectcolor', {'k'});
intersectlinestyle = ft_getopt(varargin, 'intersectlinestyle', {'-'});
intersectlinewidth = ft_getopt(varargin, 'intersectlinewidth', 2);
ncirc              = ft_getopt(varargin, 'ncirc', 15);
scalealpha         = ft_getopt(varargin, 'scalealpha', 'no');

% mesh related inputs
facecolor          = ft_getopt(varargin, 'facecolor', [0.781 0.762 0.664]);
edgecolor          = ft_getopt(varargin, 'edgecolor', 'none');
facealpha          = ft_getopt(varargin, 'facealpha', 1);
edgealpha          = ft_getopt(varargin, 'edgealpha', 0);

if rmin < 1 * ft_scalingfactor('mm', unit)
  error('cfg.rmin must be equal or larger than 1 mm');
end

if ~isempty(meshplot)
  % Mesh should be a cell-array
  if isstruct(meshplot)
    tmp = meshplot;
    meshplot = cell(size(tmp));
    for i=1:numel(tmp)
      meshplot{i} = tmp(i);
    end
  elseif iscell(meshplot)
    % do nothing
  else
    meshplot = {};
  end
  
  % replace pnt by pos
  for k = 1:numel(meshplot)
    meshplot{k} = fixpos(meshplot{k});
  end
  
  for k = 1:numel(meshplot)
    if ~isfield(meshplot{k}, 'pos') || ~isfield(meshplot{k}, 'tri')
      % error('the mesh should be a structure with pos and tri');
      meshplot{k}.pos = [];
      meshplot{k}.tri = [];
    end
  end
  
  % facecolor and edgecolor should be cell array
  if ~iscell(facecolor)
    tmp = facecolor;
    if ischar(tmp)
      facecolor = {tmp};
    elseif ismatrix(tmp) && size(tmp, 2) == 3
      facecolor = cell(size(tmp,1), 1);
      for i=1:size(tmp,1)
        facecolor{i} = tmp(i, 1:3);
      end
    else
      facecolor = {};
    end
  end
  if ~iscell(edgecolor)
    tmp = edgecolor;
    if ischar(tmp)
      edgecolor = {tmp};
    elseif ismatrix(tmp) && size(tmp, 2) == 3
      edgecolor = cell(size(tmp,1), 1);
      for i=1:size(tmp,1)
        edgecolor{i} = tmp(i, 1:3);
      end
    else
      edgecolor = {};
    end
  end
  
  % make sure each mesh has plotting options specified
  if numel(meshplot) > 1
    nmesh = numel(meshplot);
    if numel(facecolor) < numel(meshplot)
      for m = numel(facecolor)+1:nmesh
        facecolor{m} = facecolor{1};
      end
    end
    if numel(edgecolor) < numel(meshplot)
      for m = numel(edgecolor)+1:nmesh
        edgecolor{m} = edgecolor{1};
      end
    end
    if numel(facealpha) < numel(meshplot)
      for m = numel(facealpha)+1:nmesh
        facealpha(m) = facealpha(1);
      end
    end
    if numel(edgealpha) < numel(meshplot)
      for m = numel(edgealpha)+1:nmesh
        edgealpha(m) = edgealpha(1);
      end
    end
  end
end

if strcmp(sli, '2d') || strcmp(sli, '3d')
  if isempty(meshplot)
    error('plotting a slice requires a mesh as input')
  else
    dointersect = 1;
  end
else
  dointersect = 0;
end

if dointersect % check intersection inputs
  % Color and linestyle should be cell-array
  if ~iscell(intersectcolor)
    tmp = intersectcolor;
    if ischar(tmp)
      intersectcolor = {tmp};
    elseif ismatrix(tmp) && size(tmp, 2) == 3
      intersectcolor = cell(size(tmp,1), 1);
      for i=1:size(tmp,1)
        intersectcolor{i} = tmp(i, 1:3);
      end
    else
      intersectcolor = {};
    end
  end
  if ischar(intersectlinestyle)
    intersectlinestyle = {intersectlinestyle};
  elseif iscell(intersectlinestyle)
    % do nothing
  else
    intersectlinestyle = {};
  end
  
  % Make sure each intersection has plotting options specified
  if numel(meshplot) > 1
    nmesh = numel(meshplot);
    if numel(intersectcolor) < numel(meshplot)
      for m = numel(intersectcolor)+1:nmesh
        intersectcolor{m} = intersectcolor{1};
      end
    end
    if numel(intersectlinewidth) < numel(meshplot)
      for m = numel(intersectlinewidth)+1:nmesh
        intersectlinewidth(m) = intersectlinewidth(1);
      end
    end
    if numel(intersectlinestyle) < numel(meshplot)
      for m = numel(intersectlinestyle)+1:nmesh
        intersectlinestyle{m} = intersectlinestyle{1};
      end
    end
  end % end intersection plotting checks
end % end dointersect checks

if dointersect
  % Set the orientation of the slice plane
  if strcmp(ori, 'x')
    oriX = 1; oriY = 0; oriZ = 0;
  elseif strcmp(ori, 'y')
    oriX = 0; oriY = 1; oriZ = 0;
  elseif strcmp(ori, 'z')
    oriX = 0; oriY = 0; oriZ = 1;
  else
    error('ori must be "x", "y" or "z"')
  end
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

if strcmp(scalerad, 'yes')
  rmax = rmin+(radius-rmin)*radscf; % maximum radius of the clouds
else
  rmax = ones(length(pos), 1)*radius; % each cloud has the same radius
end

if dointersect
  % Generate Circle Points
  angles = linspace(0,2*pi,50);
  x = cos(angles)';
  y = sin(angles)';
  slicedim = zeros(length(angles),1);
  
  if strcmp(slicepos, 'auto') % Search each slice for largest area of data
    % Find the potential limits of the interpolation
    intxmax = max(pos(:,1))+radius; intxmin = min(pos(:,1))-radius;
    intymax = max(pos(:,2))+radius; intymin = min(pos(:,2))-radius;
    intzmax = max(pos(:,3))+radius; intzmin = min(pos(:,3))-radius;
    
    % Define potential slices with data
    if oriX; potent_slices = round(intxmin):round(intxmax); end;
    if oriY; potent_slices = round(intymin):round(intymax); end;
    if oriZ; potent_slices = round(intzmin):round(intzmax); end;
    
    area = NaN(length(pos),length(potent_slices)); % preallocate matrix of electrode interpolation area for each slice
    for s = 1:length(potent_slices) % only search slices that could potentially contain data
      distance = NaN(length(pos),1); % preallocate vector for each electrodes distance from the slice
      for c = 1:length(pos) % cloud loop
        indpos = pos(c, :);
        if oriX; distance(c) = abs(indpos(1)-potent_slices(s)); end
        if oriY; distance(c) = abs(indpos(2)-potent_slices(s)); end
        if oriZ; distance(c) = abs(indpos(3)-potent_slices(s)); end
        
        if distance(c) < rmax(c) % if there is any data from this electrode in this slice
          % find the circle points for the interpolation of this electrode
          xmax = rmax(c)*x;
          ymax = rmax(c)*y;
          
          % find the maximum radius of a cross section of the virtual sphere in the given slice
          xmaxdif = abs(xmax-distance(c));
          imindif = find(xmaxdif == min(xmaxdif), 1); % index of x value closest to distance(e)
          rcmax = abs(ymax(imindif));
          
          area(c, s) = 0.5*pi*rcmax^2;
        else % area must be zero
          area(c, s) = 0;
        end
      end % end electrode loop
    end % end slice loop
    totalarea = sum(area);
    
    slicepos = zeros(nslices,1);
    for n = 1:nslices
      imaxslice = find(totalarea == max(totalarea), 1);       % index of the slice with the maximum area
      slicepos(n) = potent_slices(imaxslice);                 % position of the yet unlisted slice with the maximum area
      totalarea(imaxslice-minspace:imaxslice+minspace) = 0;   % change the totalarea of the chosen slice and those within minspace to 0 so that it is not chosen again
    end
  end
  
  % Pre-allocate logical array specifying whether intersection with a given
  % mesh (k) actually exists within a given slice (s)
  intersect_exists = zeros(numel(slicepos), numel(meshplot));
end

% draw figure
if strcmp(sli, '2d')
  % Pre-allocate interpolation limits of each slice to facilitate
  % finding overall limits of all slices after plotting
  xsmax = NaN(numel(slicepos),1); xsmin = NaN(numel(slicepos),1);
  ysmax = NaN(numel(slicepos),1); ysmin = NaN(numel(slicepos),1);
  zsmax = NaN(numel(slicepos),1); zsmin = NaN(numel(slicepos),1);
  
  for s = 1:numel(slicepos) % slice loop
    subplot(numel(slicepos),1,s); hold on;
    
    % Pre-allocate interpolation limits of each cloud to facilitate
    % finding slice limits after plotting
    xcmax = NaN(length(pos(:,1)),1); xcmin = NaN(length(pos(:,1)),1);
    ycmax = NaN(length(pos(:,1)),1); ycmin = NaN(length(pos(:,1)),1);
    zcmax = NaN(length(pos(:,1)),1); zcmin = NaN(length(pos(:,1)),1);
    
    for c = 1:length(pos(:,1)) % cloud loop
      indpos = pos(c, :);
      % Calculate distance from slice
      if oriX; distance = abs(indpos(1)-slicepos(s)); end
      if oriY; distance = abs(indpos(2)-slicepos(s)); end
      if oriZ; distance = abs(indpos(3)-slicepos(s)); end
      
      if distance < rmax(c)
        if strcmp(slicetype, 'surf')
          xmax = rmax(c)*x;
          ymax = rmax(c)*y;
          
          if strcmp(scalealpha, 'yes')
            maxalpha = (rmax(c)-distance)/rmax(c);
          else
            maxalpha = 1;
          end
          
          % find the maximum radius of a cross section of the virtual sphere in the given slice
          xmaxdif = abs(xmax-distance);
          imindif = find(xmaxdif == min(xmaxdif), 1); % index of x value closest to distance(e)
          rcmax = abs(ymax(imindif));
          
          % Determine points along outermost circle
          xe = rcmax*x;
          ye = rcmax*y;
          
          % Jitter values of points in the slice plane so no surfaces overlap
          slicedime = slicedim+(0.01*rand*ones(length(x), 1));
          
          % Plot concentric circles
          for n = 0:ncirc-1 % circle loop
            xo = xe*((ncirc-n)/ncirc);    % outer x points
            yo = ye*((ncirc-n)/ncirc);    % outer z points
            xi = xe*((ncirc-1-n)/ncirc);  % inner x points
            yi = ye*((ncirc-1-n)/ncirc);  % inner z points
            if n == ncirc-1 % for the last concentric circle
              if oriX; hs = fill3(slicepos(s)+slicedime, indpos(2)+xo, indpos(3)+yo, val(c)); end
              if oriY; hs = fill3(indpos(1)+xo, slicepos(s)+slicedime, indpos(3)+yo, val(c)); end
              if oriZ; hs = fill3(indpos(1)+xo, indpos(2)+yo, slicepos(s)+slicedime, val(c)); end
            else
              if oriX; hs = fill3([slicepos(s)+slicedime; slicepos(s)+slicedim], [indpos(2)+xo; indpos(2)+xi], [indpos(3)+yo; indpos(3)+yi], val(c)); end
              if oriY; hs = fill3([indpos(1)+xo; indpos(1)+xi], [slicepos(s)+slicedime; slicepos(s)+slicedim], [indpos(3)+yo; indpos(3)+yi], val(c)); end
              if oriZ; hs = fill3([indpos(1)+xo; indpos(1)+xi], [indpos(2)+yo; indpos(2)+yi], [slicepos(s)+slicedime; slicepos(s)+slicedim], val(c)); end
            end
            set(hs, 'EdgeColor', 'none', 'FaceAlpha', maxalpha*n/ncirc)
          end % end circle loop
          
          % find the limits of the plotted surfaces for this electrode
          if oriX
            xcmax(c) = max(slicedime+slicepos(s)); xcmin(c) = min(slicedime+slicepos(s));
            ycmax(c) = max(xe+indpos(2)); ycmin(c) = min(xe+indpos(2));
            zcmax(c) = max(ye+indpos(3)); zcmin(c) = min(ye+indpos(3));
          elseif oriY
            xcmax(c) = max(xe+indpos(1)); xcmin(c) = min(xe+indpos(1));
            ycmax(c) = max(slicedime+slicepos(s)); ycmin(c) = min(slicedime+slicepos(s));
            zcmax(c) = max(ye+indpos(3)); zcmin(c) = min(ye+indpos(3));
          elseif oriZ
            xcmax(c) = max(xe+indpos(1)); xcmin(c) = min(xe+indpos(1));
            ycmax(c) = max(ye+indpos(2)); ycmin(c) = min(ye+indpos(2));
            zcmax(c) = max(slicedime+slicepos(s)); zcmin(c) = min(slicedime+indpos(s));
          end
        elseif strcmp(slicetype, 'point')
          rng(0, 'twister');                          % random number generator
          npoints = round(ptdens*pi*rmax(c)^2);       % number of points based on area of cloud cross section
          azimuth = 2*pi*rand(npoints,1);             % azimuthal angle for each point
          radii = rmax(c)*(rand(npoints,1).^ptgrad);  % radius value for each point
          radii = sort(radii);                        % sort radii in ascending order so they are plotted from inside out
          % convert to Carthesian; note that second input controls third output
          if oriX; [y,z,x] = sph2cart(azimuth, zeros(npoints,1)+0.01*rand(npoints,1), radii); end
          if oriY; [x,z,y] = sph2cart(azimuth, zeros(npoints,1)+0.01*rand(npoints,1), radii); end
          if oriZ; [x,y,z] = sph2cart(azimuth, zeros(npoints,1)+0.01*rand(npoints,1), radii); end
          
          % color axis with radius scaling
          if strcmp(cgrad, 'white') % color runs up to white
            fcolidx = ceil(cmid) + sign(colscf(c))*floor(abs(colscf(c)*cmid));
            if fcolidx == 0; fcolidx = 1; end
            fcol = cmapsc(fcolidx,:); % color [Nx3]
            ptcol = [linspace(fcol(1), 1, npoints)' linspace(fcol(2), 1, npoints)' linspace(fcol(3), 1, npoints)'];
          elseif isscalar(cgrad) % color runs down towards colorbar middle
            rnorm = radii/rmax(c); % normalized radius
            if radscf(c)>=.5 % extreme values
              ptcol = val(c) - (flip(1-rnorm).^inv(cgrad))*val(c); % scaled fun [Nx1]
            elseif radscf(c)<.5 % values closest to zero
              ptcol = val(c) + (flip(1-rnorm).^inv(cgrad))*abs(val(c)); % scaled fun [Nx1]
            end
          else
            error('cfg.colorgrad should be either ''white'' or a scalar determining color falloff')
          end
          
          % draw the points
          if oriX; scatter3(x+slicepos(s), y+indpos(2), z+indpos(3), ptsize, ptcol, '.'); end
          if oriY; scatter3(x+indpos(1), y+slicepos(s), z+indpos(3), ptsize, ptcol, '.'); end
          if oriZ; scatter3(x+indpos(1), y+indpos(2), z+slicepos(s), ptsize, ptcol, '.'); end
          
          % find the limits of the plotted points for this electrode
          if oriX
            xcmax(c) = max(x+slicepos(s)); xcmin(c) = min(x+slicepos(s));
            ycmax(c) = max(y+indpos(2)); ycmin(c) = min(y+indpos(2));
            zcmax(c) = max(z+indpos(3)); zcmin(c) = min(z+indpos(3));
          elseif oriY
            xcmax(c) = max(x+indpos(1)); xcmin(c) = min(x+indpos(1));
            ycmax(c) = max(y+slicepos(s)); ycmin(c) = min(y+slicepos(s));
            zcmax(c) = max(z+indpos(3)); zcmin(c) = min(z+indpos(3));
          elseif oriZ
            xcmax(c) = max(x+indpos(1)); xcmin(c) = min(x+indpos(1));
            ycmax(c) = max(y+indpos(2)); ycmin(c) = min(y+indpos(2));
            zcmax(c) = max(z+slicepos(s)); zcmin(c) = min(z+slicepos(s));
          end
          
        end % cloudtype
      end % if distance < rmax(c)
    end % cloud loop
    if dointersect
      if oriX; ori = [1 0 0]; loc = [slicepos(s) 0 0]; end
      if oriY; ori = [0 1 0]; loc = [0 slicepos(s) 0]; end
      if oriZ; ori = [0 0 1]; loc = [0 0 slicepos(s)]; end
      
      % normalise the orientation vector to one
      ori = ori./sqrt(sum(ori.^2));
      
      % shift the location to be along the orientation vector
      loc = ori*dot(loc,ori);
      
      % determine three points on the plane
      inplane = eye(3) - (eye(3) * ori') * ori;
      v1 = loc + inplane(1,:);
      v2 = loc + inplane(2,:);
      v3 = loc + inplane(3,:);
      
      for k = 1:numel(meshplot)
        
        % only plot if the mesh actually intersects the plane
        xmmax = max(meshplot{k}.pos(:,1)); xmmin = min(meshplot{k}.pos(:,1));
        ymmax = max(meshplot{k}.pos(:,2)); ymmin = min(meshplot{k}.pos(:,2));
        zmmax = max(meshplot{k}.pos(:,3)); zmmin = min(meshplot{k}.pos(:,3));
        
        if oriX
          if slicepos(s) < xmmax && slicepos(s) > xmmin
            intersect_exists(s,k) = 1;
          else
            intersect_exists(s,k) = 0;
          end
        elseif oriY
          if slicepos(s) < ymmax && slicepos(s) > ymmin
            intersect_exists(s,k) = 1;
          else
            intersect_exists(s,k) = 0;
          end
        elseif oriZ
          if slicepos(s) < zmmax && slicepos(s) > zmmin
            intersect_exists(s,k) = 1;
          else
            intersect_exists(s,k) = 0;
          end
        end
        
        if intersect_exists(s,k)
          [xmesh, ymesh, zmesh] = intersect_plane(meshplot{k}.pos, meshplot{k}.tri, v1, v2, v3);
          
          % draw each individual line segment of the intersection
          if ~isempty(xmesh)
            p = patch(xmesh', ymesh', zmesh', nan(1, size(xmesh,1)));
            if ~isempty(intersectcolor),     set(p, 'EdgeColor', intersectcolor{k}); end
            if ~isempty(intersectlinewidth), set(p, 'LineWidth', intersectlinewidth(k)); end
            if ~isempty(intersectlinestyle), set(p, 'LineStyle', intersectlinestyle{k}); end
          end
          
          % find the limits of the lines and add them to the limits of the
          % interpolation to facilitate finding the limits of the slice
          xcmax(end+1) = max(xmesh(:)); xcmin(end+1) = min(xmesh(:));
          ycmax(end+1) = max(ymesh(:)); ycmin(end+1) = min(ymesh(:));
          zcmax(end+1) = max(zmesh(:)); zcmin(end+1) = min(zmesh(:));
        end
      end % end mesh loop
    end % end if dointersect
    
    % Find limits of this particular slice
    xsmax(s) = max(xcmax); xsmin(s) = min(xcmin);
    ysmax(s) = max(ycmax); ysmin(s) = min(ycmin);
    zsmax(s) = max(zcmax); zsmin(s) = min(zcmin);
    
    % Color Settings
    colormap(cmap);
    if ~isempty(clim) && clim(2)>clim(1)
      caxis(gca, clim);
    end
    
    % Axis and View Settings
    set(gca, 'DataAspectRatio', [1 1 1])
    if oriX
      view([90 0]);
    elseif oriY
      view([180 0]);
    elseif oriZ
      view([90 90]);
    end
    
    % Add Title to Differentiate Slices
    if oriX; title(['slicepos = [' num2str(slicepos(s)) ' 0 0]']); end
    if oriY; title(['slicepos = [0 ' num2str(slicepos(s)) ' 0]']); end
    if oriZ; title(['slicepos = [0 0 ' num2str(slicepos(s)) ']']); end
  end
  
  % Set matching limits in the non-slice dimensions for each slice
  for s = 1:numel(slicepos) % slice loop
    subplot(numel(slicepos),1,s);
    if oriX
      xlim([xsmin(s)-2 xsmax(s)+2]);
      ylim([min(ysmin)-2 max(ysmax)+2]);
      zlim([min(zsmin)-2 max(zsmax)+2]);
    elseif oriY
      xlim([min(xsmin)-2 max(xsmax)+2]);
      ylim([ysmin(s)-2 ysmax(s)+2]);
      zlim([min(zsmin)-2 max(zsmax)+2]);
    elseif oriZ
      xlim([min(xsmin)-2 max(xsmax)+2]);
      ylim([min(ysmin)-2 max(ysmax)+2]);
      zlim([zsmin(s)-2 zsmax(s)+2]);
    end
  end
  
else % plot 3d cloud
  % generate point cloud(s)
  hold on;
  for n = 1:size(pos,1) % cloud loop
    % point cloud with radius scaling
    rng(0, 'twister'); % random number generator
    npoints   = round(ptdens*(4/3)*pi*rmax(n)^3);     % number of points based on cloud volume
    elevation = asin(2*rand(npoints,1)-1);            % elevation angle for each point
    azimuth   = 2*pi*rand(npoints,1);                 % azimuth angle for each point
    radii     = rmax(n)*(rand(npoints,1).^ptgrad);    % radius value for each point
    radii     = sort(radii);                          % sort radii in ascending order so they are plotted from inside out
    [x,y,z]   = sph2cart(azimuth, elevation, radii);  % convert to Carthesian
    
    % color axis with radius scaling
    if strcmp(cgrad, 'white')                   % color runs up to white
      indx  = ceil(cmid) + sign(colscf(n))*floor(abs(colscf(n)*cmid));
      indx  = max(min(indx,size(cmapsc,1)),1);  % index should fall within the colormap
      fcol  = cmapsc(indx,:);                   % color [Nx3]
      ptcol = [linspace(fcol(1), 1, npoints)' linspace(fcol(2), 1, npoints)' linspace(fcol(3), 1, npoints)'];
    elseif isscalar(cgrad)                      % color runs down towards colorbar middle
      rnorm = radii/rmax(n);                    % normalized radius
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
  
  if ~isempty(meshplot)
    for k = 1:numel(meshplot) % mesh loop
      ft_plot_mesh(meshplot{k}, 'facecolor', facecolor{k}, 'EdgeColor', edgecolor{k}, ...
        'facealpha', facealpha(k), 'edgealpha', edgealpha(k));
      material dull
    end % end mesh loop
    
    if dointersect % plot the outlines on the mesh
      for s = 1:numel(slicepos) % slice loop
        if oriX; ori = [1 0 0]; loc = [slicepos(s) 0 0]; end
        if oriY; ori = [0 1 0]; loc = [0 slicepos(s) 0]; end
        if oriZ; ori = [0 0 1]; loc = [0 0 slicepos(s)]; end
        
        % normalise the orientation vector to one
        ori = ori./sqrt(sum(ori.^2));
        
        % shift the location to be along the orientation vector
        loc = ori*dot(loc,ori);
        
        % determine three points on the plane
        inplane = eye(3) - (eye(3) * ori') * ori;
        v1 = loc + inplane(1,:);
        v2 = loc + inplane(2,:);
        v3 = loc + inplane(3,:);
        
        for k = 1:numel(meshplot)
          
          % only plot if the mesh actually intersects the plane
          xmmax = max(meshplot{k}.pos(:,1)); xmmin = min(meshplot{k}.pos(:,1));
          ymmax = max(meshplot{k}.pos(:,2)); ymmin = min(meshplot{k}.pos(:,2));
          zmmax = max(meshplot{k}.pos(:,3)); zmmin = min(meshplot{k}.pos(:,3));
          
          if oriX
            if slicepos(s) < xmmax && slicepos(s) > xmmin
              intersect_exists(s,k) = 1;
            else
              intersect_exists(s,k) = 0;
            end
          elseif oriY
            if slicepos(s) < ymmax && slicepos(s) > ymmin
              intersect_exists(s,k) = 1;
            else
              intersect_exists(s,k) = 0;
            end
          elseif oriZ
            if slicepos(s) < zmmax && slicepos(s) > zmmin
              intersect_exists(s,k) = 1;
            else
              intersect_exists(s,k) = 0;
            end
          end
          
          if intersect_exists(s,k)
            [xmesh, ymesh, zmesh] = intersect_plane(meshplot{k}.pos, meshplot{k}.tri, v1, v2, v3);
            
            % draw each individual line segment of the intersection
            if ~isempty(xmesh)
              p = patch(xmesh', ymesh', zmesh', nan(1, size(xmesh,1)));
              if ~isempty(intersectcolor),     set(p, 'EdgeColor', intersectcolor{k}); end
              if ~isempty(intersectlinewidth), set(p, 'LineWidth', intersectlinewidth(k)); end
              if ~isempty(intersectlinestyle), set(p, 'LineStyle', intersectlinestyle{k}); end
            end
          end
        end % end mesh loop
      end % end slice loop
    end % end plotting outline
  end % end mesh plotting
  
  % axis settings
  axis off
  axis vis3d
  axis equal
  
  % Color settings
  colormap(cmap);
  if ~isempty(clim) && clim(2)>clim(1)
    caxis(gca, clim);
  end
end
