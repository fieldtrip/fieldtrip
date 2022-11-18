function hs = ft_plot_headshape(headshape, varargin)

% FT_PLOT_HEADSHAPE visualizes the shape of a head from a variety of
% acquisition system. Usually the head shape is measured with a
% Polhemus tracker and some proprietary software (e.g. from CTF, BTi
% or Yokogawa). The headshape and fiducials can be used for coregistration.
% If present in the headshape, the location of the fiducials will also
% be shown.
%
% Use as
%   ft_plot_headshape(shape, ...)
% where the shape is a structure obtained from FT_READ_HEADSHAPE.
%
% Optional arguments should come in key-value pairs and can include
%   'vertexcolor'  = color specification as [r g b] values or a string, for example 'brain', 'cortex', 'skin', 'red', 'r'
%   'vertexsize'   = scalar value specifying the size of the vertices (default = 10)
%   'fidcolor'     = color specification as [r g b] values or a string, for example 'brain', 'cortex', 'skin', 'red', 'r'
%   'fidmarker'    = ['.', '*', '+',  ...]
%   'fidlabel'     = ['yes', 'no', 1, 0, 'true', 'false']
%   'transform'    = transformation matrix for the fiducials, converts MRI voxels into head shape coordinates
%   'unit'         = string, convert to the specified geometrical units (default = [])
%   'axes'         = boolean, whether to plot the axes of the 3D coordinate system (default = false)
%
% Example
%   shape = ft_read_headshape(filename);
%   ft_plot_headshape(shape)
%
% See also FT_PLOT_MESH, FT_PLOT_ORTHO

% Copyright (C) 2009, Cristiano Micheli
% Copyright (C) 2009-2022, Robert Oostenveld
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

% rename pnt into pos
headshape = fixpos(headshape);

if ~isstruct(headshape) && isnumeric(headshape) && size(headshape,2)==3
  % the input seems like a list of points, convert into something that resembles a headshape
  ws1 = warning('off', 'MATLAB:warn_r14_stucture_assignment');
  headshape.pos = headshape;
  warning(ws1);
end

% the default behavior depends on whether there is a triangulated surface or not
hastri = isfield(headshape, 'tri');

% get the optional input arguments
if hastri
  vertexcolor = ft_getopt(varargin, 'vertexcolor',  'none');
  facecolor   = ft_getopt(varargin, 'facecolor',    [1 1 1]/2);
  edgecolor   = ft_getopt(varargin, 'edgecolor',    'none');
else
  vertexcolor = ft_getopt(varargin, 'vertexcolor',  'r');
  facecolor   = ft_getopt(varargin, 'facecolor',    'none');
  edgecolor   = ft_getopt(varargin, 'edgecolor',    'none');
end
vertexsize    = ft_getopt(varargin, 'vertexsize',   10);
material_     = ft_getopt(varargin, 'material');            % do not confuse with /Applications/MATLAB_R2020b.app/toolbox/matlab/graph3d/material.m
tag           = ft_getopt(varargin, 'tag',         '');
fidcolor      = ft_getopt(varargin, 'fidcolor',     'g');
fidmarker     = ft_getopt(varargin, 'fidmarker',    '*');
fidlabel      = ft_getopt(varargin, 'fidlabel',     true);
transform     = ft_getopt(varargin, 'transform');
unit          = ft_getopt(varargin, 'unit');
axes_         = ft_getopt(varargin, 'axes', false);         % do not confuse with built-in (/Applications/MATLAB_R2020b.app/toolbox/matlab/graphics/axis/axes)

if ~isempty(unit)
  headshape = ft_convert_units(headshape, unit);
end

% color management, the other colors are handled in ft_plot_mesh
if ischar(fidcolor) && exist([fidcolor '.m'], 'file')
  fidcolor = eval(fidcolor);
end

% start with empty return values
hs = [];

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

mesh = keepfields(headshape, {'pos', 'tri', 'tet', 'hex', 'color', 'unit', 'coordsys'});
h  = ft_plot_mesh(mesh, 'vertexcolor', vertexcolor, 'vertexsize', vertexsize, 'facecolor', facecolor, 'edgecolor', edgecolor, 'material', material_, 'axes', axes_, 'tag', tag);
hs = [hs; h];

if isfield(headshape, 'fid')
  fid = headshape.fid;
  if ~isempty(transform)
    % spatially transform the fiducials
    % FIXME what is the reason for this?
    fid.pos = ft_warp_apply(transform, fid.pos);
  end
  
  % show the fiducial labels
  for i=1:size(fid.pos,1)
    h  = plot3(fid.pos(i,1), fid.pos(i,2), fid.pos(i,3), 'Marker', fidmarker, 'MarkerEdgeColor', fidcolor);
    hs = [hs; h];
    if isfield(fid, 'label') && istrue(fidlabel)
      % the text command does not like int or single position values
      x = double(fid.pos(i, 1));
      y = double(fid.pos(i, 2));
      z = double(fid.pos(i, 3));
      str = sprintf('%s', fid.label{i});
      h   = text(x, y, z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
      hs  = [hs; h];
    end
  end
end

if isfield(headshape, 'coordsys')
  % add a context sensitive menu to change the 3d viewpoint to top|bottom|left|right|front|back
  menu_viewpoint(gca, headshape.coordsys)
end

if nargout==0
  clear hs
end

if ~holdflag
  hold off
end

