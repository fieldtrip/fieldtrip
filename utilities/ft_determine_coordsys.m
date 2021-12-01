function [data] = ft_determine_coordsys(data, varargin)

% FT_DETERMINE_COORDSYS plots a geometrical object, allowing you to perform
% a visual check on the coordinatesystem, the units and on the anatomical
% labels for the coordinate system axes.
%
% Use as
%   [dataout] = ft_determine_coordsys(datain, ...)
% where the input data structure can be either
%  - an anatomical MRI
%  - a cortical or head surface mesh
%  - an electrode, gradiometer or optode definition
%  - a volume conduction model of the head
% or most other FieldTrip structures that represent geometrical information.
%
% Additional optional input arguments should be specified as key-value pairs
% and can include
%   interactive  = string, 'yes' or 'no' (default = 'yes')
%   axisscale    = scaling factor for the reference axes and sphere (default = 1)
%   clim         = lower and upper anatomical MRI limits (default = [0 1])
%
% This function wil pop up a figure that allows you to check whether the
% alignment of the object relative to the coordinate system axes is correct
% and what the anatomical labels of the coordinate system axes are. You
% should switch on the 3D rotation option in the figure panel to rotate and
% see the figure from all angles. To change the anatomical labels of the
% coordinate system, you should press the corresponding keyboard button.
%
% Recognized and supported coordinate systems are 'ctf', 'bti', '4d', 'yokogawa',
% 'eeglab', 'neuromag', 'itab', 'acpc', 'spm', 'mni', 'fsaverage', 'tal', 'scanras',
% 'scanlps', 'dicom'.
% 
% Furthermore, supported coordinate systems that do not specify the origin are 'ras',
% 'als', 'lps', etc. See https://www.fieldtriptoolbox.org/faq/coordsys for more
% details.
%
% See also FT_CONVERT_COORDSYS, FT_DETERMINE_UNITS, FT_CONVERT_UNITS, FT_PLOT_AXES, FT_PLOT_XXX

% Copyright (C) 2015-2021, Jan-Mathijs Schoffelen
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

dointeractive = ft_getopt(varargin, 'interactive', 'yes');
axisscale     = ft_getopt(varargin, 'axisscale', 1);  % this is used to scale the axmax and rbol
clim          = ft_getopt(varargin, 'clim', [0 1]);   % this is used to scale the orthoplot
fontsize      = ft_getopt(varargin, 'fontsize');      % this is passed to ft_plot_axes

data  = ft_checkdata(data, 'hasunit', 'yes');
dtype = ft_datatype(data);

% the high-level data structures are detected with ft_datatype, but there are
% also some low-level data structures that need to be supproted here
if strcmp(dtype, 'unknown')
  if isfield(data, 'fid') || (isfield(data, 'tri') && isfield(data, 'pos'))
    dtype = 'headshape';
  elseif isfield(data, 'hex') && isfield(data, 'pos')
    dtype = 'mesh';
  elseif isfield(data, 'tet') && isfield(data, 'pos')
    dtype = 'mesh';
  elseif ~strcmp(ft_headmodeltype(data), 'unknown')
    dtype = 'headmodel';
  elseif ~strcmp(ft_senstype(data), 'unknown')
    dtype = 'sens';
  end
elseif strcmp(dtype, 'mesh+label')
  % we don't care about the labels here
  dtype = 'mesh';
end

if isfield(data, 'coordsys') && ~isempty(data.coordsys)
  % ensure that it is in lower case
  data.coordsys = lower(data.coordsys);
  % print the interpretation of the coordinate system
  [labelx, labely, labelz] = coordsys2label(data.coordsys, 2, 0);
  fprintf('The positive x-axis is pointing towards %s\n', labelx);
  fprintf('The positive y-axis is pointing towards %s\n', labely);
  fprintf('The positive z-axis is pointing towards %s\n', labelz);
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
      ft_error('don''t know which volumetric parameter to plot');
    end

    % the volumetric data needs to be interpolated onto three orthogonal planes
    % determine a resolution that is close to, or identical to the original resolution
    [corner_vox, corner_head] = cornerpoints(data.dim, data.transform);
    diagonal_head = norm(range(corner_head));
    diagonal_vox  = norm(range(corner_vox));
    resolution    = (diagonal_head+eps)/(diagonal_vox+eps); % this is in units of "data.unit"
    
    % scale funparam between 0 and 1
    if ~isa(funparam, 'double') % avoid integer datatypes to allow for scaling
      funparam = double(funparam);
    end
    dmin = min(funparam(:));
    dmax = max(funparam(:));
    funparam  = (funparam-dmin)./(dmax-dmin);
    
    clear ft_plot_slice
    ft_plot_ortho(funparam, 'transform', data.transform, 'unit', data.unit, 'resolution', resolution, 'style', 'intersect', 'clim', clim);
    axis vis3d
    view([110 36]);
    
  case 'source'
    if isfield(data, 'inside') && ~isfield(data, 'tri')
      % only plot the source locations that are inside the volume conduction model
      ft_plot_mesh(data.pos(data.inside, :));
    else
      ft_plot_mesh(data, 'edgecolor','none', 'facecolor', [0.6 0.8 0.6], 'facealpha', 0.6);
    end
    camlight;

  case 'dip'
    ft_plot_mesh(data, 'edgecolor','none', 'facecolor', 'none');
    camlight;

  case 'headshape'
    ft_plot_headshape(data);
    camlight;

  case {'mesh', 'source+mesh'}
    ft_plot_mesh(data);
    camlight;

  case 'headmodel'
    ft_plot_headmodel(data);
    camlight;

  case {'grad' 'elec' 'opto' 'sens'}
    ft_plot_sens(data, 'label', 'label');
    camlight;

  case {'raw', 'timelock', 'freq', 'mvar', 'freqmvar', 'comp'}
    % the data may contain a gradiometer or electrode definition
    if isfield(data, 'grad')
      ft_plot_sens(data.grad);
    elseif isfield(data, 'elec')
      ft_plot_sens(data.elec, 'label', 'label');
    elseif isfield(data, 'opto')
      ft_plot_sens(data.opto, 'label', 'label');
    end

  case 'unknown'
end % switch dtype{k}

if isfield(data, 'tri')
  % this makes the 3-D object easier to understand
  camlight
  lighting gouraud
end

% plot the 3-D axes, labels, and sphere at the origin
ft_plot_axes(data, 'axisscale', axisscale, 'fontsize', fontsize);

if istrue(dointeractive)

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
  while ~any(strcmp(origin, {'a', 'i', 's', 'n'}))
    origin = input('Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), s(scanner origin), n(ot a landmark)? ', 's');
  end

  % some coordinate systems are identical or very similar, see https://www.fieldtriptoolbox.org/faq/coordsys
  if origin=='a' && strcmp(orientation, 'ras')
    coordsys = 'acpc'; % also used for spm, mni, tal
  elseif origin=='i' && strcmp(orientation, 'als')
    coordsys = 'ctf'; % also used for 4d, bti, eeglab
  elseif origin=='i' && strcmp(orientation, 'ras')
    coordsys = 'neuromag'; % also used for itab
  elseif origin=='s' && strcmp(orientation, 'ras')
    coordsys = 'scanras'; % also used for nifti
  elseif origin=='s' && strcmp(orientation, 'lps')
    coordsys = 'scanlps'; % also used for dicom
  else
    % only use the orientation, not the origin
    coordsys = orientation;
  end

  data.coordsys = coordsys;
end % if interactive
