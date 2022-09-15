function hs = ft_plot_sens(sens, varargin)

% FT_PLOT_SENS visualizes the EEG, MEG or NIRS sensor array.
%
% Use as
%   ft_plot_sens(sens, ...)
% where the first argument is the sensor array as returned by FT_READ_SENS or
% by FT_PREPARE_VOL_SENS.
%
% Optional input arguments should come in key-value pairs and can include
%   'label'           = show the label, can be 'off', 'label', 'number' (default = 'off')
%   'chantype'        = string or cell-array with strings, for example 'meg' (default = 'all')
%   'unit'            = string, convert the sensor array to the specified geometrical units (default = [])
%   'axes'            = boolean, whether to plot the axes of the 3D coordinate system (default = false)
%   'fontcolor'       = string, color specification (default = 'k')
%   'fontsize'        = number, sets the size of the text (default = 10)
%   'fontunits'       =
%   'fontname'        =
%   'fontweight'      =
%
% The following options apply to MEG magnetometers and/or gradiometers
%   'coil'            = true/false, plot each individual coil (default = false)
%   'orientation'     = true/false, plot a line for the orientation of each coil (default = false)
%   'coilshape'       = 'point', 'circle', 'square', 'sphere', or 'disc' (default is automatic)
%   'coilsize'        = diameter or edge length of the coils (default is automatic)
% The following options apply to EEG electrodes
%   'elec'            = true/false, plot each individual electrode (default = false)
%   'orientation'     = true/false, plot a line for the orientation of each electrode (default = false)
%   'elecshape'       = 'point', 'circle', 'square', 'sphere', or 'disc' (default is automatic)
%   'elecsize'        = diameter of the electrodes (default is automatic)
%   'headshape'       = headshape, required for elecshape 'disc'
% The following options apply to NIRS optodes
%   'opto'            = true/false, plot each individual optode (default = false)
%   'orientation'     = true/false, plot a line for the orientation of each optode (default = false)
%   'optoshape'       = 'point', 'circle', 'square', 'sphere', or 'disc' (default is automatic)
%   'optosize'        = diameter of the optodes (default is automatic)
%
% The following options apply when electrodes/coils/optodes are NOT plotted individually
%   'style'           = plotting style for the points representing the channels, see plot3 (default = [])
%   'marker'          = marker type representing the channels, see plot3 (default = '.')
% The following options apply when electrodes/coils/optodes are plotted individually
%   'facecolor'       = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', or an Nx3 or Nx1 array where N is the number of faces (default is automatic)
%   'edgecolor'       = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', color of channels or coils (default is automatic)
%   'facealpha'       = transparency, between 0 and 1 (default = 1)
%   'edgealpha'       = transparency, between 0 and 1 (default = 1)
%
% Example
%   sens = ft_read_sens('Subject01.ds', 'senstype', 'meg');
%   figure; ft_plot_sens(sens, 'coilshape', 'point', 'style', 'r*')
%   figure; ft_plot_sens(sens, 'coilshape', 'circle')
%   figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', true, 'chantype', 'meggrad')
%   figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', false, 'orientation', true)
%
% See also FT_READ_SENS, FT_PLOT_HEADSHAPE, FT_PLOT_HEADMODEL

% Copyright (C) 2009-2022, Robert Oostenveld, Arjen Stolk
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

% ensure that the sensor description is up-to-date
sens = ft_datatype_sens(sens);

% get the optional input arguments
label           = ft_getopt(varargin, 'label', 'off');
chantype        = ft_getopt(varargin, 'chantype');
unit            = ft_getopt(varargin, 'unit');
axes_           = ft_getopt(varargin, 'axes', false);     % do not confuse with built-in function
orientation     = ft_getopt(varargin, 'orientation', false);
% these have to do with the font
fontcolor       = ft_getopt(varargin, 'fontcolor', 'k');  % default is black
fontsize        = ft_getopt(varargin, 'fontsize',   get(0, 'defaulttextfontsize'));
fontname        = ft_getopt(varargin, 'fontname',   get(0, 'defaulttextfontname'));
fontweight      = ft_getopt(varargin, 'fontweight', get(0, 'defaulttextfontweight'));
fontunits       = ft_getopt(varargin, 'fontunits',  get(0, 'defaulttextfontunits'));

% this is for MEG magnetometer and/or gradiometer arrays
coil            = ft_getopt(varargin, 'coil', false);
coilshape       = ft_getopt(varargin, 'coilshape'); % default depends on the input, see below
coilsize        = ft_getopt(varargin, 'coilsize');  % default depends on the input, see below
% this is for EEG electrode arrays
elec            = ft_getopt(varargin, 'elec', false);
elecshape       = ft_getopt(varargin, 'elecshape'); % default depends on the input, see below
elecsize        = ft_getopt(varargin, 'elecsize');  % default depends on the input, see below
headshape       = ft_getopt(varargin, 'headshape', []); % needed for elecshape 'disc'
% this is for NIRS optode arrays
opto            = ft_getopt(varargin, 'opto', false);
optoshape       = ft_getopt(varargin, 'optoshape'); % default depends on the input, see below
optosize        = ft_getopt(varargin, 'optosize');  % default depends on the input, see below

iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');
isnirs = ft_senstype(sens, 'nirs');

% make sure that the options are consistent with the data
if iseeg
  individual = elec;
  sensshape  = elecshape;
  senssize   = elecsize;
elseif ismeg
  individual = coil;
  sensshape  = coilshape;
  senssize   = coilsize;
elseif isnirs
  % this has not been tested
  individual = opto;
  sensshape  = optoshape;
  senssize   = optosize;
else
  ft_warning('unknown sensor array description');
  individual = false;
  sensshape  = [];
  senssize   = [];
end

% this is simply passed to plot3
style           = ft_getopt(varargin, 'style');
marker          = ft_getopt(varargin, 'marker', '.');

% this is simply passed to ft_plot_mesh
if strcmp(sensshape, 'sphere') || strcmp(sensshape, 'disc')
  edgecolor     = ft_getopt(varargin, 'edgecolor', 'none');
else
  edgecolor     = ft_getopt(varargin, 'edgecolor', 'k');
end
facecolor       = ft_getopt(varargin, 'facecolor');  % default depends on the input, see below
facealpha       = ft_getopt(varargin, 'facealpha',   1);
edgealpha       = ft_getopt(varargin, 'edgealpha',   1);

if ischar(chantype)
  % this should be a cell-array
  chantype = {chantype};
end

if ~isempty(ft_getopt(varargin, 'coilorientation'))
  % for backward compatibility, added on 17 Aug 2016
  ft_warning('the coilorientation option is deprecated, please use "orientation" instead')
  orientation = ft_getopt(varargin, 'coilorientation');
end

if ~isempty(ft_getopt(varargin, 'coildiameter'))
  % for backward compatibility, added on 6 July 2016
  % the senssize is the diameter for a circle, or the edge length for a square
  ft_warning('the coildiameter option is deprecated, please use "coilsize" instead')
  senssize = ft_getopt(varargin, 'coildiameter');
end

if ~isempty(unit)
  % convert the sensor description to the specified units
  sens = ft_convert_units(sens, unit);
end

if isempty(sensshape)
  if ft_senstype(sens, 'neuromag')
    if strcmp(chantype, 'megmag')
      sensshape = 'point'; % these cannot be plotted as squares
    else
      sensshape = 'square';
    end
  elseif ft_senstype(sens, 'meg')
    sensshape = 'circle';
  else
    sensshape = 'point';
  end
end

if isempty(senssize)
  % start with a size expressed in millimeters
  switch ft_senstype(sens)
    case 'neuromag306'
      senssize = 30; % FIXME this is only an estimate
    case 'neuromag122'
      senssize = 35; % FIXME this is only an estimate
    case 'ctf151'
      senssize = 20;
    case 'ctf275'
      senssize = 18;
    otherwise
      if strcmp(sensshape, 'sphere') || strcmp(sensshape, 'disc')
        senssize = 4; % assuming spheres/discs are used for intracranial electrodes, diameter is about 4mm
      elseif strcmp(sensshape, 'point')
        senssize = 30;
      else
        senssize = 10;
      end
  end
  % convert from mm to the units of the sensor array
  senssize = senssize/ft_scalingfactor(sens.unit, 'mm');
end

% color management
if isempty(facecolor) % set default color depending on shape
  if strcmp(sensshape, 'point')
    facecolor = 'k';
  elseif strcmp(sensshape, 'circle') || strcmp(sensshape, 'square')
    facecolor = 'none';
  elseif strcmp(sensshape, 'sphere') || strcmp(sensshape, 'disc')
    facecolor = 'b';
  end
end
if ischar(facecolor) && exist([facecolor '.m'], 'file')
  facecolor = feval(facecolor);
end
if ischar(edgecolor) && exist([edgecolor '.m'], 'file')
  edgecolor = feval(edgecolor);
end

% select a subset of channels and coils to be plotted
if ~isempty(chantype)
  % remove the balancing from the sensor definition, e.g. 3rd order gradients, PCA-cleaned data or ICA projections
  sens = undobalancing(sens);
  
  chansel = match_str(sens.chantype, chantype);
  
  % remove the channels that are not selected
  sens.label    = sens.label(chansel);
  sens.chanpos  = sens.chanpos(chansel,:);
  sens.chantype = sens.chantype(chansel);
  sens.chanunit = sens.chanunit(chansel);
  if isfield(sens, 'chanori')
    % this is only present for MEG sensor descriptions
    sens.chanori  = sens.chanori(chansel,:);
  end
  
  % remove the magnetometer and gradiometer coils that are not in one of the selected channels
  if isfield(sens, 'tra') && isfield(sens, 'coilpos')
    sens.tra     = sens.tra(chansel,:);
    coilsel      = any(sens.tra~=0,1);
    sens.coilpos = sens.coilpos(coilsel,:);
    sens.coilori = sens.coilori(coilsel,:);
    sens.tra     = sens.tra(:,coilsel);
  end
  
  % FIXME note that I have not tested this on any complicated electrode definitions
  % remove the electrodes that are not in one of the selected channels
  if isfield(sens, 'tra') && isfield(sens, 'elecpos')
    sens.tra     = sens.tra(chansel,:);
    elecsel      = any(sens.tra~=0,1);
    sens.elecpos = sens.elecpos(elecsel,:);
    sens.tra     = sens.tra(:,elecsel);
  end
  
end % selecting channels and coils

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

if istrue(individual)
  % get the position of all individual coils, electrodes or optodes
  if isfield(sens, 'coilpos')
    pos = sens.coilpos;
  elseif isfield(sens, 'elecpos')
    pos = sens.elecpos;
  elseif isfield(sens, 'optopos')
    pos = sens.optopos;
  end
  if isfield(sens, 'coilori')
    ori = sens.coilori;
  elseif isfield(sens, 'elecori')
    ori = sens.elecori;
  elseif isfield(sens, 'optoori')
    ori = sens.optoori;
  else
    ori = [];
  end
  
else
  % determine the position of each channel, which is for example the mean of
  % two bipolar electrodes, or the bottom coil of a axial gradiometer, or
  % the center between two coils of a planar gradiometer
  if isfield(sens, 'chanpos')
    pos = sens.chanpos;
  else
    pos = [];
  end
  if isfield(sens, 'chanori')
    ori = sens.chanori;
  else
    ori = [];
  end
  
end % if istrue(individual)

if isempty(ori)
  if ~isempty(headshape)
    % the following code uses some functions from the computer vision toolbox
    % ft_hastoolbox('vision', -1);
    
    % how many local points on the headshape are used for estimating the local norm
    npoints = 25;
    
    % calculate local norm vectors
    for i=1:size(pos,1)
      % compute the distance to all headshape points
      d = sqrt( (pos(i,1)-headshape.pos(:,1)).^2 + (pos(i,2)-headshape.pos(:,2)).^2 + (pos(i,3)-headshape.pos(:,3)).^2 );
      [dum, idx] = sort(d);
      x = headshape.pos(idx(1:npoints),1);
      y = headshape.pos(idx(1:npoints),2);
      z = headshape.pos(idx(1:npoints),3);
      ptCloud = pointCloud([x y z]);
      nrm = pcnormals(ptCloud);
      u = nrm(:,1);
      v = nrm(:,2);
      w = nrm(:,3);
      
      % compute the headshape center
      C = mean(headshape.pos,1);
      % the vector should be pointing away from the center, otherwise flip it
      for k = 1:numel(x)
        p1 = C - [x(k) y(k) z(k)];
        p2 =     [u(k) v(k) w(k)];
        angle = atan2(norm(cross(p1,p2)),p1*p2');
        if ~(angle > pi/2 || angle < -pi/2)
          u(k) = -u(k);
          v(k) = -v(k);
          w(k) = -w(k);
        end
      end
      Fn = nanmean([u v w],1);
      Fn = Fn * (1/sqrt(sum(Fn.^2,2))); % normalize
      ori(i,:) = Fn;
    end % for
    
  elseif ~any(isnan(pos(:))) && size(pos,1)>2
    % determine orientations based on a surface triangulation of the sensors
    tri = projecttri(pos, 'delaunay');
    ori = normals(pos, tri);
    
  elseif size(pos,1)>4
    % determine orientations by fitting a sphere to the sensors
    try
      tmp = pos(~any(isnan(pos), 2),:); % remove rows that contain a nan
      center = fitsphere(tmp);
    catch
      center = [nan nan nan];
    end
    for i=1:size(pos,1)
      ori(i,:) = pos(i,:) - center;
      ori(i,:) = ori(i,:)/norm(ori(i,:));
    end
  else
    ori = nan(size(pos));
  end
  
end % if empty(ori)

if mean(isnan(ori(:)))>0.25
  % more than a quarter of the sensor orientations cannot be determined
  % the ones that have been determined probably don't make sense either
  if iseeg
    ft_notice('orienting all EEG electrodes along the z-axis')
  elseif ismeg
    ft_notice('orienting all MEG sensors along the z-axis')
  elseif isnirs
    ft_notice('orienting all NIRS optodes along the z-axis')
  end
  ori(:,1) = 0;
  ori(:,2) = 0;
  ori(:,3) = 1;
elseif any(isnan(ori(:)))
  % only some of the sensor positions cannot be determined
  % the others probably do make sense
  if iseeg
    ft_notice('orienting some EEG electrodes along the z-axis')
  elseif ismeg
    ft_notice('orienting some MEG sensors along the z-axis')
  elseif isnirs
    ft_notice('orienting some NIRS optodes along the z-axis')
  end
  ori(isnan(ori(:,1)),1) = 0;
  ori(isnan(ori(:,2)),2) = 0;
  ori(isnan(ori(:,3)),3) = 1;
end

if istrue(orientation)
  scale = ft_scalingfactor('mm', sens.unit)*20; % draw a line segment of 20 mm
  for i=1:size(pos,1)
    x = [pos(i,1) pos(i,1)+ori(i,1)*scale];
    y = [pos(i,2) pos(i,2)+ori(i,2)*scale];
    z = [pos(i,3) pos(i,3)+ori(i,3)*scale];
    line(x, y, z)
  end
end

switch sensshape
  case 'point'
    if ~isempty(style)
      % the style can include the color and/or the shape of the marker
      % check whether the marker shape is specified
      possible = {'+', 'o', '*', '.', 'x', 'v', '^', '>', '<', 'square', 'diamond', 'pentagram', 'hexagram'};
      specified = false(size(possible));
      for i=1:numel(possible)
        specified(i) = ~isempty(strfind(style, possible{i}));
      end
      if any(specified)
        % the marker shape is specified in the style option
        hs = plot3(pos(:,1), pos(:,2), pos(:,3), style, 'MarkerSize', senssize);
      else
        % the marker shape is not specified in the style option, use the marker option instead and assume that the style option represents the color
        hs = plot3(pos(:,1), pos(:,2), pos(:,3), 'Marker', marker, 'MarkerSize', senssize, 'Color', style, 'Linestyle', 'none');
      end
    else
      % the style is not specified, use facecolor for the marker
      % if the marker is '.' it will show points that do not depend on the size, in all other cases (e.g. 'o') the size is relevant
      hs = scatter3(pos(:,1), pos(:,2), pos(:,3), senssize.^2, facecolor, marker);
    end
    
  case 'sphere'
    plotsens(pos, ori, [], senssize, sensshape, 'edgecolor', edgecolor, 'facecolor', facecolor, 'edgealpha', edgealpha, 'facealpha', facealpha);
    
  case 'disc'
    plotsens(pos, ori, [], senssize, sensshape, 'edgecolor', edgecolor, 'facecolor', facecolor, 'edgealpha', edgealpha, 'facealpha', facealpha);
    
  case 'circle'
    plotsens(pos, ori, [], senssize, sensshape, 'edgecolor', edgecolor, 'facecolor', facecolor, 'edgealpha', edgealpha, 'facealpha', facealpha);
    
  case 'square'
    % determine the rotation-around-the-axis of each sensor
    % this is only applicable for neuromag planar gradiometers
    if ft_senstype(sens, 'neuromag')
      [nchan, ncoil] = size(sens.tra);
      chandir = nan(nchan,3);
      for i=1:nchan
        poscoil = find(sens.tra(i,:)>0);
        negcoil = find(sens.tra(i,:)<0);
        if numel(poscoil)==1 && numel(negcoil)==1
          % planar gradiometer
          direction = sens.coilpos(poscoil,:)-sens.coilpos(negcoil,:);
          direction = direction/norm(direction);
          chandir(i,:) = direction;
        elseif (numel([poscoil negcoil]))==1
          % magnetometer
        elseif numel(poscoil)>1 || numel(negcoil)>1
          ft_error('cannot work with balanced gradiometer definition')
        end
      end
    else
      chandir = [];
    end
    
    plotsens(pos, ori, chandir, senssize, sensshape, 'edgecolor', edgecolor, 'facecolor', facecolor, 'edgealpha', edgealpha, 'facealpha', facealpha);
    
  otherwise
    ft_error('incorrect shape');
end % switch

if ~isempty(label) && ~any(strcmp(label, {'off', 'no'}))
  
  % determine the amount of offset for the labels
  if strcmp(sensshape, 'point')
    % determine the median of the distance to the nearest neighbour
    sensdist = triu(dist(sens.chanpos'),1);
    sensdist(sensdist==0) = Inf;
    sensdist = min(sensdist,[], 2);
    sensdist = median(sensdist);
    % the offset is based on distance between sensors
    offset = 0.5 * sensdist;
  else
    % the offset is based on size of the sensors
    offset = 1.5 * senssize;
  end
  
  if isinf(offset)
    % this happens in case there is only one sensor and the size has not been specified
    offset = ft_scalingfactor('mm', sens.unit)*10; % displace the label by 10 mm
  end
  
  for i=1:size(pos,1)
    switch label
      case {'on', 'yes', 'label', 'labels'}
        if ~individual
          str = sens.label{i};
        elseif ismeg
          % individual MEG coils never have a label
          str = '';
        elseif iseeg
          if isequal(sens.chanpos, sens.elecpos)
            % the names of the electrodes and channels can be interchanged
            str = sens.label{i};
          else
            % the names of the individual electrodes are not known
            str = '';
          end
        elseif isnirs
          % optodes have individual names
          str = sens.optolabel{i};
        end
      case {'number' 'numbers'}
        str = num2str(i);
      otherwise
        ft_error('unsupported value for option ''label''');
    end % switch
    % shift the label with a certain offset
    x = pos(i,1) + offset * ori(i,1);
    y = pos(i,2) + offset * ori(i,2);
    z = pos(i,3) + offset * ori(i,3);
    text(x, y, z, str, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight, 'horizontalalignment', 'center', 'verticalalignment', 'middle');
  end % for each channel
end % if label

axis vis3d
axis equal

if istrue(axes_)
  % plot the 3D axes, this depends on the units and coordsys
  ft_plot_axes(sens);
end

if isfield(sens, 'coordsys')
  % add a context sensitive menu to change the 3d viewpoint to top|bottom|left|right|front|back
  menu_viewpoint(gca, sens.coordsys)
end

if ~holdflag
  hold off
end

if ~nargout
  clear hs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION all optional inputs are passed to ft_plot_mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'%%%%%%%%%%%%%%%%%%
function plotsens(senspos, sensori, sensdir, senssize, sensshape, varargin)

% start with a single template coil at [0 0 0], oriented towards [0 0 1]
switch sensshape
  case 'sphere'
    [pos, tri] = mesh_sphere(100);
    pos(:,1) = pos(:,1)/2;  % unit diameter
    pos(:,2) = pos(:,2)/2;  % unit diameter
    pos(:,3) = pos(:,3)/2;  % unit diameter
  case 'disc'
    [pos, tri] = mesh_cylinder(36, 2);
    pos(:,1) = pos(:,1)/2;  % unit diameter
    pos(:,2) = pos(:,2)/2;  % unit diameter
    pos(:,3) = pos(:,3)/10;
  case 'circle'
    pos = circle(24);
    tri = [];
  case 'square'
    pos = square;
    tri = [];
end

nsens = size(senspos,1);
npos  = size(pos,1);

mesh.pos  = nan(nsens*npos,3);
mesh.poly = nan(nsens, npos);   % this will be used for the edge of the coil or square
mesh.tri  = nan(0, 3);          % this will be used for the discs

% determine the scaling of the coil as homogenous transformation matrix
s = scale([senssize senssize senssize]);

for i=1:nsens
  x  = sensori(i,1);
  y  = sensori(i,2);
  z  = sensori(i,3);
  ph = atan2(y, x)*180/pi;
  th = atan2(sqrt(x^2+y^2), z)*180/pi;
  % determine the rotation and translation of the coil as homogenous transformation matrix
  r1 = rotate([0 th 0]);
  r2 = rotate([0 0 ph]);
  t  = translate(senspos(i,:));
  
  % determine the initial rotation of the coil as homogenous transformation matrix
  if isempty(sensdir)
    % none of the coils needs to be rotated around their axis, this applies to circular coils
    r0 = eye(4);
  elseif ~all(isfinite(sensdir(i,:)))
    % the rotation around the axis of this coil is not known
    r0 = nan(4);
  else
    % express the direction of sensitivity of the planar channel relative to the orientation of the channel
    dir = ft_warp_apply(inv(r2*r1), sensdir(i,:));
    x = dir(1);
    y = dir(2);
    % determine the rotation
    rh = atan2(y, x)*180/pi;
    r0 = rotate([0 0 rh]);
  end
  
  switch sensshape
    case {'sphere' 'disc'}
      % construct a single mesh with separate triangles for all sensors
      sel = ((i-1)*npos+1):(i*npos);
      mesh.pos(sel,:) = ft_warp_apply(t*r2*r1*r0*s, pos);
      mesh.tri        = cat(1, mesh.tri, tri + (i-1)*npos);
      mesh.poly       = [];
    case {'circle' 'square'}
      % construct a single mesh with separate polygons for all sensors
      sel = ((i-1)*npos+1):(i*npos);
      mesh.pos(sel,:) = ft_warp_apply(t*r2*r1*r0*s, pos); % scale, rotate and translate the template coil vertices, skip the central vertex
      mesh.poly(i,:)  = sel;                              % this is a polygon connecting all edge points
      mesh.tri        = [];
  end
  
end % for each sensor

% use either poly or tri for plotting
if isempty(mesh.tri)
  mesh = rmfield(mesh, 'tri');
end
if isempty(mesh.poly)
  mesh = rmfield(mesh, 'poly');
end

% plot all polygons together
ft_plot_mesh(mesh, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION return a circle with unit diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos] = circle(n)
phi = linspace(0, 2*pi, n+1)';
x = cos(phi);
y = sin(phi);
z = zeros(size(phi));
pos = [x y z]/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION return a square with unit edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos] = square
pos = [
   0.5  0.5 0
  -0.5  0.5 0
  -0.5 -0.5 0
   0.5 -0.5 0
   0.5  0.5 0 % this closes the square
  ];
