function ft_plot_axes(object, varargin)

% FT_PLOT_AXES adds three axes of 150 mm and a 10 mm sphere at the origin to the
% present 3-D figure. The axes and sphere are scaled according to the units of the
% geometrical object that is passed to this function. Furthermore, when possible,
% the axes labels will represent the aanatomical labels corresponding to the
% specified coordinate system.
%
% Use as
%   ft_plot_axes(object)
%
% Additional optional input arguments should be specified as key-value pairs
% and can include
%   'axisscale'    = scaling factor for the reference axes and sphere (default = 1)
%   'unit'         = string, convert the data to the specified geometrical units (default = [])
%   'coordsys'     = string, assume the data to be in the specified coordinate system (default = 'unknown')
%   'fontcolor'    = string, color specification (default = [1 .5 0], i.e. orange)
%   'fontsize'     = number, sets the size of the text (default is automatic)
%   'fontunits'    =
%   'fontname'     =
%   'fontweight'   =
%
% See also FT_PLOT_SENS, FT_PLOT_MESH, FT_PLOT_ORTHO, FT_PLOT_HEADSHAPE, FT_PLOT_DIPOLE, FT_PLOT_VOL

% Copyright (C) 2015, Jan-Mathijs Schoffelen
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

axisscale = ft_getopt(varargin, 'axisscale', 1);   % this is used to scale the axmax and rbol
unit      = ft_getopt(varargin, 'unit');
coordsys  = ft_getopt(varargin, 'coordsys');
% these have to do with the font
fontcolor   = ft_getopt(varargin, 'fontcolor', [1 .5 0]); % default is orange
fontsize    = ft_getopt(varargin, 'fontsize',   get(0, 'defaulttextfontsize'));
fontname    = ft_getopt(varargin, 'fontname',   get(0, 'defaulttextfontname'));
fontweight  = ft_getopt(varargin, 'fontweight', get(0, 'defaulttextfontweight'));
fontunits   = ft_getopt(varargin, 'fontunits',  get(0, 'defaulttextfontunits'));

% color management
if ischar(fontcolor) && exist([fontcolor '.m'], 'file')
  fontcolor = eval(fontcolor);
end

if ~isempty(object) && ~isempty(unit)
  % convert the object to the specified units
  object = ft_convert_units(object, unit);
elseif ~isempty(object) &&  isempty(unit)
  % take the units from the object
  object = ft_determine_units(object);
  unit = object.unit;
elseif  isempty(object) && ~isempty(unit)
  % there is no object, but the units have been specified
elseif  isempty(object) &&  isempty(unit)
  ft_warning('units are not known, not plotting axes')
  return
end

if ~isempty(object) && ~isfield(object, 'coordsys')
  % set it to unknown, that makes the subsequent code easier
  object.coordsys = 'unknown';
end

if ~isempty(object) && ~isempty(coordsys)
  % check the user specified coordinate system with the one in the object
  assert(strcmp(coordsys, unit.coordsys), 'coordsys is inconsistent with the object')
elseif ~isempty(object) &&  isempty(coordsys)
  % take the coordinate system from the object
  coordsys = object.coordsys;
elseif  isempty(object) && ~isempty(coordsys)
  % there is no object, but the coordsys has been specified
elseif  isempty(object) &&  isempty(coordsys)
  % this is not a problem per see
  coordsys = 'unknown';
end

axmax = 150 * ft_scalingfactor('mm', unit);
rbol  =   5 * ft_scalingfactor('mm', unit);

% this is useful if the anatomy is from a non-human primate or rodent
axmax = axisscale*axmax;
rbol  = axisscale*rbol;

fprintf('The axes are %g %s long in each direction\n', axmax, unit);
fprintf('The diameter of the sphere at the origin is %g %s\n', 2*rbol, unit);

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

prevhold = ishold;
hold on

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
[O.pos, O.tri] = icosahedron42;
O.pos = O.pos.*rbol;
ft_plot_mesh(O, 'edgecolor', 'none');

% create the labels that are to be plotted along the axes
[labelx, labely, labelz] = coordsys2label(coordsys, 3, 1);

% add the labels to the axis
text(xdat(1,1), ydat(1,1), zdat(1,1), labelx{1}, 'linewidth', 2, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
text(xdat(1,2), ydat(1,2), zdat(1,2), labely{1}, 'linewidth', 2, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
text(xdat(1,3), ydat(1,3), zdat(1,3), labelz{1}, 'linewidth', 2, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
text(xdat(2,1), ydat(2,1), zdat(2,1), labelx{2}, 'linewidth', 2, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
text(xdat(2,2), ydat(2,2), zdat(2,2), labely{2}, 'linewidth', 2, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);
text(xdat(2,3), ydat(2,3), zdat(2,3), labelz{2}, 'linewidth', 2, 'color', fontcolor, 'fontunits', fontunits, 'fontsize', fontsize, 'fontname', fontname, 'fontweight', fontweight);

if ~prevhold
  hold off
end
