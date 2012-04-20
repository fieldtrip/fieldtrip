function [data] = ft_determine_coordsys(data, varargin)

% FT_DETERMINE_COORDSYS plots a geometrical object, allowing you to perform
% a visual check on the coordinatesystem, the units and on the anatomical
% labels for the coordinate system axes.
%
% Use as
%   [dataout] = ft_determine_coordsys(datain, 'key1', value1, ...)
% where the input data structure can be
%  - an anatomical MRI
%  - an electrode or gradiometer definition
%  - a volume conduction model
% or most other FieldTrip structures that represent geometrical information.
%
% Additional optional input arguments should be specified as key-value pairs 
% and can include
%   interactive  = string, 'yes' or 'no' (default = 'yes')
%
% This function wil pop up a figure that allows you to check whether the
% alignment of the object relative to the coordinate system axes is correct
% and what the anatomical labels of the coordinate system axes are. You should
% switch on the 3D rotation option in the figure panel to rotate and see the
% figure from all angles. To change the anatomical labels of the coordinate
% system, you should press the corresponding keyboard button.
%
% See also FT_VOLUMEREALIGN, FT_VOLUMERESLICE

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

dointeractive = ft_getopt(varargin, 'interactive', 'yes');
axisscale     = ft_getopt(varargin, 'axisscale', 1); % this is used to scale the axmax and rbol

data  = ft_checkdata(data);
dtype = ft_datatype(data);
data  = ft_convert_units(data);
unit  = data.unit;

% the high-level data structures are detected with ft_datatype, but there are
% also some low-level data structures that need to be supproted here
if strcmp(dtype, 'unknown')
  if isfield(data, 'fid') || (isfield(data, 'tri') && isfield(data, 'pnt'))
    dtype = 'headshape';
  elseif ~strcmp(ft_voltype(data), 'unknown')
    dtype = 'headmodel';
  elseif ~strcmp(ft_senstype(data), 'unknown')
    dtype = 'sens';
  end
end

% determine the size of the "unit" sphere in the origin and the length of the axes
switch unit
  case 'mm'
    axmax = 150;
    rbol  = 5;
  case 'cm'
    axmax = 15;
    rbol  = 0.5;
  case 'm'
    axmax = 0.15;
    rbol  = 0.005;
  otherwise
    error('unknown units (%s)', unit);
end

% this is useful if the anatomy is from a non-human primate or rodent
axmax = axisscale*axmax;
rbol  = axisscale*rbol;

if isfield(data, 'coordsys') && ~isempty(data.coordsys)
  label = cell(3,1);
  if length(data.coordsys)==3 && length(intersect(data.coordsys, 'rlasif'))==3
    for i=1:3
      switch data.coordsys(i)
        case 'l'
          label{i} = 'the left';
        case 'r'
          label{i} = 'the right';
        case 'i'
          label{i} = 'inferior';
        case 's'
          label{i} = 'superior';
        case 'a'
          label{i} = 'anterior';
        case 'p'
          label{i} = 'posterior';
        otherwise
          error('incorrect letter in the coordsys');
      end % switch
    end % for each of the three axes
  elseif strcmpi(data.coordsys, 'itab') || strcmpi(data.coordsys, 'neuromag')
    label{1} = 'the right';
    label{2} = 'anterior';
    label{3} = 'superior';
  elseif strcmpi(data.coordsys, 'ctf') || strcmpi(data.coordsys, '4d') || strcmpi(data.coordsys, 'bti')
    label{1} = 'anterior';
    label{2} = 'the left';
    label{3} = 'superior';
  elseif strcmpi(data.coordsys, 'tal') || strcmpi(data.coordsys, 'mni') || strcmpi(data.coordsys, 'spm')
    label{1} = 'the right';
    label{2} = 'anterior';
    label{3} = 'superior';
  else
    error('unsupported coordsys');
  end
  
  fprintf('The positive x-axis is pointing towards %s\n', label{1});
  fprintf('The positive y-axis is pointing towards %s\n', label{2});
  fprintf('The positive z-axis is pointing towards %s\n', label{3});
end

% create the labels that are to be plotted along the axes
if isfield(data, 'coordsys')
  [labelx, labely, labelz] = xyz2label(data.coordsys);
else
  [labelx, labely, labelz] = xyz2label('unknown');
end

% plot the geometrical object
% the plotting style depends on the data content
figure;
switch dtype
  case 'volume'
    funparam = [];
    if isfield(data, 'anatomy')
      funparam = data.anatomy;
    elseif isfield(data, 'gray')
      funparam = data.gray;
    elseif isfield(data, 'white')
      funparam = data.white;
    elseif isfield(data, 'brick0')
      funparam = data.brick0; % used for an atlas
    elseif isfield(data, 'brick1')
      funparam = data.brick1; % used for an atlas
    else
      % try to determine it automatically
      fn = fieldnames(data);
      for i=1:length(fn)
        if isequal(size(data.(fn{i})), data.dim)
          funparam = data.(fn{i});
          break;
        end
      end
    end
    
    if isempty(funparam)
      error('don''t know which volumetric parameter to plot');
    end
    
    % the volumetric data needs to be interpolated onto three orthogonal planes
    % determine a resolution that is close to, or identical to the original resolution
    [corner_vox, corner_head] = cornerpoints(data.dim, data.transform);
    diagonal_head = norm(range(corner_head));
    diagonal_vox  = norm(range(corner_vox));
    resolution    = diagonal_head/diagonal_vox; % this is in units of "data.unit"
    
    clear ft_plot_slice
    ft_plot_ortho(funparam, 'transform', data.transform, 'resolution', resolution, 'style', 'intersect');
    axis vis3d
    view([110 36]);
    
  case 'source'
    ft_plot_mesh(data, 'edgecolor','none', 'facecolor', [0.6 0.8 0.6], 'facealpha', 0.6);
    camlight;
    
  case 'dip'
    ft_plot_mesh(data, 'edgecolor','none', 'facecolor', 'none');
    camlight;
    
  case 'headshape'
    ft_plot_headshape(data);
    camlight;
    
  case 'headmodel'
    ft_plot_vol(data);
    camlight;
    
  case 'sens'
    ft_plot_sens(data);
    camlight;
    
  case {'raw', 'timelock', 'freq', 'mvar', 'freqmvar', 'comp'}
    % the data may contain a gradiometer or electrode definition
    if isfield(data, 'grad')
      ft_plot_sens(data.grad);
    elseif isfield(data, 'elec')
      ft_plot_sens(data.elec);
    end
    
  case 'unknown'
end % switch dtype{k}

% get the xyz-axes
xdat  = [-axmax 0 0; axmax 0 0];
ydat  = [0 -axmax 0; 0 axmax 0];
zdat  = [0 0 -axmax; 0 0 axmax];

% get the xyz-axes dotted
xdatdot = (-axmax:(axmax/15):axmax);
xdatdot = xdatdot(1:floor(numel(xdatdot)/2)*2);
xdatdot = reshape(xdatdot, [2 numel(xdatdot)/2]);
n       = size(xdatdot,2);
ydatdot = [zeros(2,n) xdatdot zeros(2,n)];
zdatdot = [zeros(2,2*n) xdatdot];
xdatdot = [xdatdot zeros(2,2*n)];

% plot axes
hl = line(xdat, ydat, zdat);
set(hl(1), 'linewidth', 1, 'color', 'r');
set(hl(2), 'linewidth', 1, 'color', 'g');
set(hl(3), 'linewidth', 1, 'color', 'b');
hld = line(xdatdot, ydatdot, zdatdot);
for k = 1:n
  set(hld(k    ), 'linewidth', 3, 'color', 'r');
  set(hld(k+n*1), 'linewidth', 3, 'color', 'g');
  set(hld(k+n*2), 'linewidth', 3, 'color', 'b');
end

% create the ball at the origin
[O.pnt, O.tri] = icosahedron42;
O.pnt = O.pnt.*rbol;
ft_plot_mesh(O, 'edgecolor', 'none');

% add the labels to the axis
text(xdat(1,1),ydat(1,1),zdat(1,1),labelx{1},'color','y','fontsize',15,'linewidth',2);
text(xdat(1,2),ydat(1,2),zdat(1,2),labely{1},'color','y','fontsize',15,'linewidth',2);
text(xdat(1,3),ydat(1,3),zdat(1,3),labelz{1},'color','y','fontsize',15,'linewidth',2);
text(xdat(2,1),ydat(2,1),zdat(2,1),labelx{2},'color','y','fontsize',15,'linewidth',2);
text(xdat(2,2),ydat(2,2),zdat(2,2),labely{2},'color','y','fontsize',15,'linewidth',2);
text(xdat(2,3),ydat(2,3),zdat(2,3),labelz{2},'color','y','fontsize',15,'linewidth',2);

if dointeractive,
  
  if ~isfield(data, 'coordsys') || isempty(data.coordsys)
    % default is yes
    value = smartinput('Do you want to change the anatomical labels for the axes [Y, n]? ', 'y');
  else
    % default is no
    value = smartinput('Do you want to change the anatomical labels for the axes [y, N]? ', 'n');
  end
  
  if strcmpi(value, 'n')
    return
  end
  
  % interactively determine orientation
  orientation = '   ';
  while ~any(strcmp(orientation(1), {'r', 'l', 'a', 'p', 's', 'i'}))
    orientation(1) = smartinput('What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? ', '');
  end
  while ~any(strcmp(orientation(2), {'r', 'l', 'a', 'p', 's', 'i'}))
    orientation(2) = smartinput('What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? ', '');
  end
  while ~any(strcmp(orientation(3), {'r', 'l', 'a', 'p', 's', 'i'}))
    orientation(3) = smartinput('What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? ', '');
  end
  
  % interactively determine origin
  origin = ' ';
  while ~any(strcmp(origin, {'a', 'i', 'n'}))
    origin = input('Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), n(ot a landmark)? ', 's');
  end
  
  if origin=='a' && strcmp(orientation, 'ras')
    coordsys = 'spm';
  elseif origin=='i' && strcmp(orientation, 'als')
    coordsys = 'ctf';
  elseif origin=='i' && strcmp(orientation, 'ras')
    coordsys = 'neuromag'; % also used for itab
  else
    % just use the orientation
    coordsys = orientation;
  end
  
  data.coordsys = coordsys;
end % if interactive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to go from aplrsi to better interpretable format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [labelx, labely, labelz] = xyz2label(str)

if ~isempty(str) && ~strcmp(str, 'unknown')
  % the first part is important for the orientations
  % the second part optionally contains information on the origin
  strx = tokenize(str, '_');
  
  switch lower(strx{1})
    case {'ras' 'itab' 'neuromag' 'spm' 'mni'}
      labelx = {'-X (left)'      '+X (right)'   };
      labely = {'-Y (posterior)' '+Y (anterior)'};
      labelz = {'-Z (inferior)'  '+Z (superior)'};
    case {'als' 'ctf' '4d', 'bti'}
      labelx = {'-X (posterior)' '+X (anterior)'};
      labely = {'-Y (right)'     '+Y (left)'};
      labelz = {'-Z (inferior)'  '+Z (superior)'};
    otherwise
      error('unknown coordsys');
  end
  
else
  labelx = {'-X (unknown)' '+X (unknown)'};
  labely = {'-Y (unknown)' '+Y (unknown)'};
  labelz = {'-Z (unknown)' '+Z (unknown)'};
end
