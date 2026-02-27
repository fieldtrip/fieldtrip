function [cfg, sensor] = ft_sensorplacement(cfg, headshape)

% FT_SENSORPLACEMENT positions sensor sensor holders over the surface of the scalp (when
% wearing a flexible cap) or over the surface of a rigid 3D-printed helmet. It uses
% a model of the sensor sensor and sensor holder, copies this repeatedly, and positions
% and orients it for each desired sensor position. The sensor positions are automatically
% determined based on a template distribution, for example the 10-20 electrode placement
% scheme or a equidistant placement scheme, but you can also provide your own sensor
% positions.
%
% Use as
%   [cfg, sensors] = ft_sensorplacement(cfg, headshape)
% where the headshape represents the scalp surface from FT_PREPARE_MESH or
% FT_MESHREALIGN. This function returns a structure array with a number of
% meshes representing the sensor sensors or sensor holders that can be plotted.
%
% The input configuration structure can contain the following
%   cfg.template      = string, filename with the STL model of the sensor or sensor holder
%   cfg.write         = 'no' or 'yes', write the sensors to STL files
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.channel       = cell-array, selection of electrode locations at which to place an sensor sensor
%   cfg.outwardshift  = number, amount to shift the sensors outward from the surface
%   cfg.rotx          = Nx1 vector with the rotation around the x-axis (default is automatic)
%   cfg.roty          = Nx1 vector with the rotation around the y-axis (default is automatic)
%   cfg.rotz          = Nx1 vector with the rotation around the z-axis (default is automatic)
%   cfg.grad          = structure with a single OPM sensor, see FT_DATATYPE_SENS
%
% The output configuration structure contains the rotations that were performed,
% which can be adjusted and used in a second iteration. The output sensor structure
% array contains the geometrical description of all sensors, following rotation and
% translation.
%
% See also FT_ELECTRODEPLACEMENT, FT_PREPARE_MESH, FT_MESHREALIGN, FT_DEFACEMESH

% Copyright (C) 2026, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    headshape
ft_preamble provenance headshape

% ensure that the input data is valid for this function, this will also do
headshape = ft_checkdata(headshape, 'datatype', 'mesh', 'feedback', 'yes');

% get the sensor positions, read them from file if needed
elec = ft_fetch_sens(cfg);

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'elec', 'template'});

% set the defaults
cfg.channel       = ft_getopt(cfg, 'channel', 'all');
cfg.write         = ft_getopt(cfg, 'write', 'no');
cfg.outwardshift  = ft_getopt(cfg, 'outwardshift', 0);
cfg.orientation   = ft_getopt(cfg, 'orientation', 'surface');
cfg.rotx          = ft_getopt(cfg, 'rotx');
cfg.roty          = ft_getopt(cfg, 'roty');
cfg.rotz          = ft_getopt(cfg, 'rotz');

% select the desired electrode positions
elec = keepfields(elec, {'elecpos', 'elecori', 'label'});
cfg.channel = ft_channelselection(cfg.channel, elec.label);
[sel1, sel2] = match_str(cfg.channel, elec.label); % sort them according to cfg.channel
elec.label   = elec.label(sel2);
elec.elecpos = elec.elecpos(sel2,:);
if isfield(elec, 'elecori')
  elec.elecori = elec.elecori(sel2,:);
else
  % use the direction perpendicular to the headshape, see below
end
nsens = length(elec.label);

% use the same rotation for each channel, or nan if not specified
if isscalar(cfg.rotx)
  cfg.rotx = ones(nsens,1) * cfg.rotx;
elseif isempty(cfg.rotx)
  cfg.rotx = nan(nsens,1);
end
if isscalar(cfg.roty)
  cfg.roty = ones(nsens,1) * cfg.roty;
elseif isempty(cfg.roty)
  cfg.roty = nan(nsens,1);
end
if isscalar(cfg.rotz)
  cfg.rotz = ones(nsens,1) * cfg.rotz;
elseif isempty(cfg.rotz)
  cfg.rotz = nan(nsens,1);
end

% read the template STL model, assume them to be in milimeter
if ischar(cfg.template)
  template = ft_read_headshape(cfg.template, 'unit', 'unknown');
  template.unit = 'mm';
else
  % use the template object as specified, it can be a grad structure with a single OPM sensor 
  template = cfg.template;
end

% project the electrodes onto the headshape surface
[dum, elec.elecpos] = project_elec(elec.elecpos, headshape.pos, headshape.tri);

if ~isfield(elec, 'elecori')
  % compute the orientation of electrodes
  elec.elecori = normals_elec(elec.elecpos, headshape.pos, headshape.tri);
end % if not elecori

for i=1:nsens
  % first shift it away from the surface
  t1 = translate([0, 0, cfg.outwardshift]);

  % determine the required orientation and rotation
  x = elec.elecori(i,1);
  y = elec.elecori(i,2);
  z = elec.elecori(i,3);

  % then rotate about z by angle γ (yaw)
  % then rotate about y by angle β (pitch)
  % then rotate about x by angle α (roll)

  if ~isnan(cfg.rotz(i))
    gamma = cfg.rotz(i)*pi/180; % convert from degrees to radians
  else
    gamma = 0;
  end

  if ~isnan(cfg.roty(i))
    beta = cfg.roty(i)*pi/180; % convert from degrees to radians
  else
    beta = asin(x);
  end

  if ~isnan(cfg.rotx(i))
    alpha = cfg.rotx(i)*pi/180; % convert from degrees to radians
  elseif x==+1 && y==0 && z==0
    alpha = 0;
    beta  = 0;
  elseif x==-1 && y==0 && z==0
    alpha = 0;
    beta  = pi;
  else
    alpha = atan2(-y, z);
  end

  % convert from radians to degrees
  gamma = gamma*180/pi;
  beta  = beta*180/pi;
  alpha = alpha*180/pi;

  % rotate around z, then around y, then around x
  r = rotate([alpha, beta, gamma]);

  % determine the final translation towards the electrode position
  t2 = translate(elec.elecpos(i,:));

  % translate, rotate, translate once more
  sensor(i) = ft_transform_geometry(t2 * r * t1, template);

  % remember the rotations
  cfg.rotz(i) = gamma;
  cfg.roty(i) = beta;
  cfg.rotx(i) = alpha;

end % for each sensor

if istrue(cfg.write)
  for i=1:nsens
    [p, f, x] = fileparts(cfg.template);
    filename = [f '_' elec.label{i} '.stl'];
    ft_info('writing %s', filename)
    ft_write_headshape(filename, sensor(i), 'format', 'stl');
  end % for each sensor
end

ft_postamble debug
