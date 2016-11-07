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
%   axisscale    = scaling factor for the reference axes and sphere (default = 1)
%   fontcolor    = string (default = 'y')
%   'unit'       = string, convert the data to the specified geometrical units (default = [])
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
fontcolor = ft_getopt(varargin, 'fontcolor', 'y'); % default is yellow
unit      = ft_getopt(varargin, 'unit');

if ~isempty(unit)
  % convert the sensor description to the specified units
  object = ft_convert_units(object, unit);
end

if isfield(object, 'unit')
  unit = object.unit;
else
  warning('units are not known, not plotting axes')
  return
end

if isfield(object, 'coordsys')
  coordsys = object.coordsys;
else
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
[labelx, labely, labelz] = xyz2label(coordsys);

% add the labels to the axis
text(xdat(1,1), ydat(1,1), zdat(1,1), labelx{1}, 'color', fontcolor, 'fontsize', 15, 'linewidth', 2);
text(xdat(1,2), ydat(1,2), zdat(1,2), labely{1}, 'color', fontcolor, 'fontsize', 15, 'linewidth', 2);
text(xdat(1,3), ydat(1,3), zdat(1,3), labelz{1}, 'color', fontcolor, 'fontsize', 15, 'linewidth', 2);
text(xdat(2,1), ydat(2,1), zdat(2,1), labelx{2}, 'color', fontcolor, 'fontsize', 15, 'linewidth', 2);
text(xdat(2,2), ydat(2,2), zdat(2,2), labely{2}, 'color', fontcolor, 'fontsize', 15, 'linewidth', 2);
text(xdat(2,3), ydat(2,3), zdat(2,3), labelz{2}, 'color', fontcolor, 'fontsize', 15, 'linewidth', 2);

if ~prevhold
  hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%
% NOTE this should be kept consistent with the longer axes labels in FT_DETERMINE_COORDSYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [labelx, labely, labelz] = xyz2label(str)

if ~isempty(str) && ~strcmp(str, 'unknown')
  % the first part is important for the orientations
  % the second part optionally contains information on the origin
  strx = tokenize(str, '_');
  
  switch lower(strx{1})
    case {'ras' 'itab' 'neuromag' 'spm' 'mni' 'tal'}
      labelx = {'-X (left)'      '+X (right)'   };
      labely = {'-Y (posterior)' '+Y (anterior)'};
      labelz = {'-Z (inferior)'  '+Z (superior)'};
    case {'als' 'ctf' '4d', 'bti'}
      labelx = {'-X (posterior)' '+X (anterior)'};
      labely = {'-Y (right)'     '+Y (left)'};
      labelz = {'-Z (inferior)'  '+Z (superior)'};
    case {'paxinos'}
      labelx = {'-X (left)'      '+X (right)'};
      labely = {'-Y (inferior)'  '+Y (superior)'};
      labelz = {'-Z (anterior)'  '+Z (posterior)'};
    case {'lps'}
      labelx = {'-X (right)'      '+X (left)'};
      labely = {'-Y (anterior)'  '+Y (posterior)'};
      labelz = {'-Z (inferior)'  '+Z (superior)'};
    otherwise
      error('unknown coordsys');
  end
  
else
  labelx = {'-X (unknown)' '+X (unknown)'};
  labely = {'-Y (unknown)' '+Y (unknown)'};
  labelz = {'-Z (unknown)' '+Z (unknown)'};
end
