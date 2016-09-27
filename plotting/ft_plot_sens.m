function hs = ft_plot_sens(sens, varargin)

% FT_PLOT_SENS plots the position of all channels or coils that comprise the EEG or
% MEG sensor array description
%
% Use as
%   ft_plot_sens(sens, ...)
% where the first argument is the sensor array as returned by FT_READ_SENS or
% by FT_PREPARE_VOL_SENS.
%
% Optional input arguments should come in key-value pairs and can include
%   'chantype'        = string or cell-array with strings, for example 'meg' (default = 'all')
%   'unit'            = string, convert the data to the specified geometrical units (default = [])
%   'label'           = show the label, can be 'off', 'label', 'number' (default = 'off')
%   'coil'            = true/false, plot each individual coil (default = false)
%   'orientation'     = true/false, plot a line for the orientation of each coil (default = false)
%
% The following option only applies when individual coils are being plotted
%   'coilsize'        = diameter or edge length of the coils (default is automatic)
%   'coilshape'       = 'point', 'circle' or 'square' (default is automatic)
%   'facecolor'       = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r', or an Nx1 array where N is the number of faces (default = 'none')
%   'edgecolor'       = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r' (default = 'k')
%   'facealpha'       = transparency, between 0 and 1 (default = 1)
%   'edgealpha'       = transparency, between 0 and 1 (default = 1)
%
% The following option only applies when individual coils are NOT being plotted
%   'style'           = plotting style for the points representing the channels, see plot3 (default = 'k.')
%
% Example
%   sens = ft_read_sens('Subject01.ds');
%   ft_plot_sens(sens, 'style', 'r*')
%
% See also FT_READ_SENS, FT_DATATYPE_SENS, FT_PLOT_HEADSHAPE, FT_PREPARE_VOL_SENS

% Copyright (C) 2009-2016, Robert Oostenveld
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

ws = warning('on', 'MATLAB:divideByZero');

% ensure that the sensor description is up-to-date (Aug 2011)
sens = ft_datatype_sens(sens);

% get the optional input arguments
label           = ft_getopt(varargin, 'label', 'off');
chantype        = ft_getopt(varargin, 'chantype');
unit            = ft_getopt(varargin, 'unit');
orientation     = ft_getopt(varargin, 'orientation', false);
% this is for MEG magnetometer and/or gradiometer arrays
coil            = ft_getopt(varargin, 'coil', false);
coilsize        = ft_getopt(varargin, 'coilsize');  % default depends on the input, see below
coilshape       = ft_getopt(varargin, 'coilshape'); % default depends on the input, see below
% this is simply passed to plot3
style           = ft_getopt(varargin, 'style', 'k.');
% this is simply passed to ft_plot_mesh
edgecolor       = ft_getopt(varargin, 'edgecolor', 'k');
facecolor       = ft_getopt(varargin, 'facecolor', 'none');
facealpha       = ft_getopt(varargin, 'facealpha',   1);
edgealpha       = ft_getopt(varargin, 'edgealpha',   1);

if ischar(chantype)
  % should be a cell array
  chantype = {chantype};
end

if ~isempty(ft_getopt(varargin, 'coilorientation'))
  % for backward compatibility, added on 17 Aug 2016
  warning('the coilorientation option is deprecated, please use orientation instead')
  orientation = ft_getopt(varargin, 'coilorientation');
end

if ~isempty(ft_getopt(varargin, 'coildiameter'))
  % for backward compatibility, added on 6 July 2016
  % the size is the diameter for a circle, or the edge length for a square
  warning('the coildiameter option is deprecated, please use coilsize instead')
  coilsize = ft_getopt(varargin, 'coilsize');
end

if ~isempty(unit)
  % convert the sensor description to the specified units
  sens = ft_convert_units(sens, unit);
end

if isempty(coilshape)
  if ft_senstype(sens, 'neuromag')
    if strcmp(chantype, 'megmag')
      coilshape = 'point'; % these cannot be plotted as squares
    else
      coilshape = 'square';
    end
  elseif ft_senstype(sens, 'meg')
    coilshape = 'circle';
  else
    coilshape = 'point';
  end
end

if isempty(coilsize)
  switch ft_senstype(sens)
    case 'neuromag306'
      coilsize = 30; % FIXME this is only an estimate
    case 'neuromag122'
      coilsize = 35; % FIXME this is only an estimate
    case 'ctf151'
      coilsize = 15; % FIXME this is only an estimate
    case 'ctf275'
      coilsize = 15; % FIXME this is only an estimate
    otherwise
      coilsize = 10;
  end
  % convert from mm to the units of the sensor array
  coilsize = coilsize/ft_scalingfactor(sens.unit, 'mm');
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

if istrue(orientation)
  if istrue(coil)
    if isfield(sens, 'coilori')
      pos = sens.coilpos;
      ori = sens.coilori;
    elseif isfield(sens, 'elecori')
      pos = sens.elecpos;
      ori = sens.elecori;
    else
      pos = [];
      ori = [];
    end
  else
    if isfield(sens, 'chanori')
      pos = sens.chanpos;
      ori = sens.chanori;
    else
      pos = [];
      ori = [];
    end
  end
  scale = ft_scalingfactor('mm', sens.unit)*20; % draw a line segment of 20 mm
  for i=1:size(pos,1)
    x = [pos(i,1) pos(i,1)+ori(i,1)*scale];
    y = [pos(i,2) pos(i,2)+ori(i,2)*scale];
    z = [pos(i,3) pos(i,3)+ori(i,3)*scale];
    line(x, y, z)
  end
end

% determine the rotation-around-the-axis of each sensor
% only applicable for neuromag planar gradiometers
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
      error('cannot work with balanced gradiometer definition')
    end
  end
end

if istrue(coil)
  % simply get the position of all coils or electrodes
  if isfield(sens, 'coilpos')
    pos = sens.coilpos;
  elseif isfield(sens, 'elecpos')
    pos = sens.elecpos;
  end
  if isfield(sens, 'coilori')
    ori = sens.coilori;
  elseif isfield(sens, 'elecori')
    ori = sens.elecori;
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
  
end % if istrue(coil)

switch coilshape
  case 'point'
    hs = plot3(pos(:,1), pos(:,2), pos(:,3), style, 'MarkerSize', 30, 'Color', edgecolor);
  case 'circle'
    plotcoil(pos, ori, [],      coilsize, coilshape, 'edgecolor', edgecolor, 'facecolor', facecolor, 'edgealpha', edgealpha, 'facealpha', facealpha);
  case 'square'
    plotcoil(pos, ori, chandir, coilsize, coilshape, 'edgecolor', edgecolor, 'facecolor', facecolor, 'edgealpha', edgealpha, 'facealpha', facealpha);
  otherwise
    error('incorrect coilshape');
end % switch

if ~isempty(label) && ~any(strcmp(label, {'off', 'no'}))
  for i=1:length(sens.label)
    switch label
      case {'on', 'yes'}
        str = sens.label{i};
      case {'label' 'labels'}
        str = sens.label{i};
      case {'number' 'numbers'}
        str = num2str(i);
      otherwise
        error('unsupported value for option ''label''');
    end % switch
    text(sens.chanpos(i,1), sens.chanpos(i,2), sens.chanpos(i,3), str);
  end % for
end % if empty or off/no

axis vis3d
axis equal

if ~nargout
  clear hs
end
if ~holdflag
  hold off
end

warning(ws); % revert to original state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION all optional inputs are passed to ft_plot_mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'%%%%%%%%%%%%%%%%%%
function plotcoil(coilpos, coilori, chandir, coilsize, coilshape, varargin)

% start with a single template coil at [0 0 0], oriented towards [0 0 1]
switch coilshape
  case 'circle'
    pos = circle(24);
  case 'square'
    pos = square;
end
ncoil = size(coilpos,1);
npos  = size(pos,1);
mesh.pos  = nan(ncoil*npos,3);
mesh.poly = nan(ncoil,npos);

% determine the scaling of the coil as homogenous transformation matrix
s  = scale([coilsize coilsize coilsize]);

for i=1:ncoil
  x  = coilori(i,1);
  y  = coilori(i,2);
  z  = coilori(i,3);
  ph = atan2(y, x)*180/pi;
  th = atan2(sqrt(x^2+y^2), z)*180/pi;
  % determine the rotation and translation of the coil as homogenous transformation matrix
  r1 = rotate([0 th 0]);
  r2 = rotate([0 0 ph]);
  t  = translate(coilpos(i,:));
  
  % determine the initial rotation of the coil as homogenous transformation matrix
  if isempty(chandir)
    % none of the coils needs to be rotated around their axis, this applies to circular coils
    r0 = eye(4);
  elseif ~all(isfinite(chandir(i,:)))
    % the rotation around the axis of this coil is not known
    r0 = nan(4);
  else
    % express the direction of sensitivity of the planar channel relative to the orientation of the channel
    dir = ft_warp_apply(inv(r2*r1), chandir(i,:));
    x = dir(1);
    y = dir(2);
    % determine the rotation
    rh = atan2(y, x)*180/pi;
    r0 = rotate([0 0 rh]);
  end
  
  % construct a single mesh with separate polygons for all coils
  sel = ((i-1)*npos+1):(i*npos);
  mesh.pos(sel,:) = ft_warp_apply(t*r2*r1*r0*s, pos); % scale, rotate and translate the template coil vertices, skip the central vertex
  mesh.poly(i,:)  = sel;                              % this is a polygon connecting all edge points
  
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
