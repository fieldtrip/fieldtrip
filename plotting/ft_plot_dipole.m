function h = ft_plot_dipole(pos, ori, varargin)

% FT_PLOT_DIPOLE makes a 3-D representation of a dipole using a sphere and a stick
% pointing along the dipole orientation
%
% Use as
%   ft_plot_dipole(pos, mom, ...)
% where pos and mom are the dipole mosition and moment. 
%
% Optional input arguments should be specified in key-value pairs and can include
%   'diameter' = number indicating sphere diameter (default = 'auto')
%   'length'   = number indicating length of the stick (default = 'auto')
%   'color'    = [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r' (default = 'r')
%   'unit'     = 'm', 'cm' or 'mm', used for automatic scaling (default = 'cm')
%   'scale'    = scale the dipole with the amplitude, can be 'none',  'both', 'diameter', 'length' (default = 'none')
%   'alpha'    = alpha value of the plotted dipole
%
% Example
%   ft_plot_dipole([0 0 0], [1 2 3], 'color', 'r', 'alpha', 1)

% Copyright (C) 2009, Robert Oostenveld
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

% get the optional input arguments
amplitudescale = ft_getopt(varargin, 'scale',    'none');
color          = ft_getopt(varargin, 'color',    'r'); % can also be a RGB triplet
diameter       = ft_getopt(varargin, 'diameter', 'auto');
length         = ft_getopt(varargin, 'length',   'auto');
unit           = ft_getopt(varargin, 'unit',     'cm');
alpha          = ft_getopt(varargin, 'alpha',     1);

% for backward compatibility, this can be changed into an error at the end of 2016
units = ft_getopt(varargin, 'units');
if ~isempty(units)
  warning('please use "unit" instead of "units"');
  unit = units;
  clear units
end
  
if isequal(diameter, 'auto')
  % the default is a 5 mm sphere
  switch unit
    case 'm'
      diameter = 0.005;
    case 'cm'
      diameter = 0.5;
    case 'mm'
      diameter = 5;
    otherwise
      error('unsupported unit');
  end
end

if isequal(length, 'auto')
  % the default is a 15 mm stick
  switch unit
    case 'm'
      length = 0.015;
    case 'cm'
      length = 1.5;
    case 'mm'
      length = 15;
    otherwise
      error('unsupported unit');
  end
end

% dipole position should be Nx3
if all(size(pos) == [3 1])
  pos = pos';
end

% dipole moment and orientation should be 3xN
if all(size(ori) == [1 3])
  ori = ori';
end

h = [];

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

% these are reused
[unitsphere.pos, unitsphere.tri] = icosahedron642;
[unitcylinder.pos, unitcylinder.tri] = cylinder(36, 2);

for i=1:size(pos,1)
  amplitude = norm(ori(:,i));
  ori(:,i) = ori(:,i) ./ amplitude;
  
  % scale the dipole diameter and length with its amplitude
  if strcmp(amplitudescale, 'length') || strcmp(amplitudescale, 'both')
    this_length = length*amplitude;
  else
    this_length = length;
  end
  if strcmp(amplitudescale, 'diameter') || strcmp(amplitudescale, 'both')
    this_diameter = diameter*amplitude;
  else
    this_diameter = diameter;
  end
  
  % start with a unit sphere and cylinder
  sphere  = unitsphere;
  stick   = unitcylinder;
  sphere.pos = ft_warp_apply(scale([0.5 0.5 0.5]), sphere.pos, 'homogeneous'); % the diameter should be 1
  stick.pos = ft_warp_apply(scale([0.5 0.5 0.5]), stick.pos, 'homogeneous'); % the length should be 1
  stick.pos = ft_warp_apply(translate([0 0 0.5]), stick.pos, 'homogeneous'); % it should start in the origin
  
  % scale the sphere
  sx = this_diameter;
  sy = this_diameter;
  sz = this_diameter;
  sphere.pos = ft_warp_apply(scale([sx sy sz]),     sphere.pos, 'homogeneous');
  
  % translate the sphere
  tx = pos(i,1);
  ty = pos(i,2);
  tz = pos(i,3);
  sphere.pos = ft_warp_apply(translate([tx ty tz]), sphere.pos, 'homogeneous');
  
  % scale the stick
  sx = this_diameter/3;
  sy = this_diameter/3;
  sz = this_length;
  stick.pos = ft_warp_apply(scale([sx sy sz]),     stick.pos, 'homogeneous');
  
  % first rotate the stick to point along the x-axis
  stick.pos = ft_warp_apply(rotate([0 90 0]),    stick.pos, 'homogeneous');
  % then rotate the stick in the desired direction
  [az, el] = cart2sph(ori(1,i), ori(2,i), ori(3,i));
  stick.pos = ft_warp_apply(rotate([0 -el*180/pi 0]),  stick.pos, 'homogeneous'); % rotate around y-axis
  stick.pos = ft_warp_apply(rotate([0  0 az*180/pi]),  stick.pos, 'homogeneous'); % rotate around z-axis
  
  % translate the stick
  tx = pos(i,1);
  ty = pos(i,2);
  tz = pos(i,3);
  stick.pos = ft_warp_apply(translate([tx ty tz]), stick.pos, 'homogeneous');
  
  % plot the sphere and the stick
  p1 = ft_plot_mesh(sphere, 'vertexcolor', 'none', 'edgecolor', false, 'facecolor', color, 'facealpha', alpha);
  h = cat(2, h(:)', p1(:)');
  clear p1;
  
  p2 = ft_plot_mesh(stick,  'vertexcolor', 'none', 'edgecolor', false, 'facecolor', color, 'facealpha', alpha);
  h = cat(2, h(:)', p2(:)');
  clear p2;
end % for each dipole

axis off
axis vis3d
axis equal

if ~holdflag
  hold off
end

if ~nargout
  clear h
end

warning(ws); %revert to original state
