function plot_dipole(pos, ori, varargin)

% PLOT_DIPOLE makes a 3-D representation of a dipole using a small sphere
% and a stick pointing along the dipole orientation
%
% Use as
%   plot_dipole(pos, mom, ...)
% where pos and mom are the dipole mosition and moment. Optional
% input arguments should be specified in key-value pairs and can
% include
%   'diameter'    number indicating sphere diameter (default = 'auto')
%   'length'      number indicating length of the stick (default = 'auto')
%   'color'       [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%   'units'       'm', 'cm' or 'mm', used for automatic scaling (default = 'cm')
%   'scale'       scale the dipole with the amplitude, can be 'none',  'both', 'diameter', 'length' (default = 'none')
%
% Example
%   plot_dipole([0 0 0], [1 2 3])

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: plot_dipole.m,v $
% Revision 1.5  2009/10/07 12:51:22  roboos
% fixed typo in help
%
% Revision 1.4  2009/09/23 14:55:36  crimic
% small fix
%
% Revision 1.3  2009/06/03 16:13:19  roboos
% updated documentation
%
% Revision 1.2  2009/06/03 09:53:18  roboos
% number of small changes, still work in progress
%
% Revision 1.1  2009/06/03 07:53:27  roboos
% first version
%

% get the optional input arguments
units       = keyval('units', varargin);       if isempty(units),         units = 'cm';        end
color       = keyval('color', varargin);       if isempty(color),         color = [1 0 0];     end
diameter    = keyval('diameter', varargin);    if isempty(diameter),      diameter = 'auto';   end
length      = keyval('length', varargin);      if isempty(length),        length = 'auto';     end
amplitudescale = keyval('scale', varargin);    if isempty(amplitudescale),  amplitudescale = 'none'; end

if isequal(diameter, 'auto')
  % the default is a 5 mm sphere
  switch units
    case 'm'
      diameter = 0.005;
    case 'cm'
      diameter = 0.5;
    case 'mm'
      diameter = 5;
    otherwise
      error('unsupported units');
  end
end

if isequal(length, 'auto')
  % the default is a 15 mm stick
  switch units
    case 'm'
      length = 0.015;
    case 'cm'
      length = 1.5;
    case 'mm'
      length = 15;
    otherwise
      error('unsupported units');
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

% everything is added to the current figure
holdflag = ishold;
hold on

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

  % create a unit sphere and cylinder
  [sphere.pnt, sphere.tri] = icosahedron642;
  sphere.pnt = warp_apply(scale([0.5 0.5 0.5]), sphere.pnt, 'homogeneous'); % the diameter should be 1
  [stick.pnt, stick.tri]   = cylinder(36, 2);
  stick.pnt = warp_apply(scale([0.5 0.5 0.5]), stick.pnt, 'homogeneous'); % the length should be 1
  stick.pnt = warp_apply(translate([0 0 0.5]), stick.pnt, 'homogeneous'); % it should start in the origin

  % scale the sphere
  sx = this_diameter;
  sy = this_diameter;
  sz = this_diameter;
  sphere.pnt = warp_apply(scale([sx sy sz]),     sphere.pnt, 'homogeneous');

  % translate the sphere
  tx = pos(i,1);
  ty = pos(i,2);
  tz = pos(i,3);
  sphere.pnt = warp_apply(translate([tx ty tz]), sphere.pnt, 'homogeneous');

  % scale the stick
  sx = this_diameter/3;
  sy = this_diameter/3;
  sz = this_length;
  stick.pnt = warp_apply(scale([sx sy sz]),     stick.pnt, 'homogeneous');

  % first rotate the stick to point along the x-axis
  stick.pnt = warp_apply(rotate([0 90 0]),    stick.pnt, 'homogeneous');
  % then rotate the stick in the desired direction
  [az, el] = cart2sph(ori(1,i), ori(2,i), ori(3,i));
  stick.pnt = warp_apply(rotate([0 -el*180/pi 0]),  stick.pnt, 'homogeneous'); % rotate around y-axis
  stick.pnt = warp_apply(rotate([0  0 az*180/pi]),  stick.pnt, 'homogeneous'); % rotate around z-axis

  % translate the stick
  tx = pos(i,1);
  ty = pos(i,2);
  tz = pos(i,3);
  stick.pnt = warp_apply(translate([tx ty tz]), stick.pnt, 'homogeneous');

  % plot the sphere and the stick
  plot_mesh(sphere, 'vertexcolor', 'none', 'edgecolor', false, 'facecolor', color);
  plot_mesh(stick,  'vertexcolor', 'none', 'edgecolor', false, 'facecolor', color);

end % for each dipole

axis off
axis vis3d
axis equal
camlight % fixme, this probably should be in the calling function

if ~holdflag
  hold off
end

