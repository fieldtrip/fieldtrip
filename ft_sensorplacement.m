function [sensor] = ft_sensorplacement(cfg, headshape)

% FT_SENSORPLACEMENT positions sensor sensor holders over the surface of the scalp (when
% wearing a flexible cap) or over the surface of a rigid 3D-printed helmet. It uses 
% a model of the sensor sensor and sensor holder, copies this repeatedly, and positions 
% and orients it for each desired sensor position. The sensor positions are automatically
% determined based on a template distribution, for example the 10-20 electrode placement
% scheme or a equidistant placement scheme, but you can also provide your own sensor
% positions.
%
% Use as
%   [sensors] = ft_sensorplacement(cfg, headshape)
% where the headshape represents the scalp surface from FT_PREPARE_MESH or
% FT_MESHREALIGN. This function returns a structure array with a number of 
% meshes representing the sensor sensors or sensor holders that can be plotted.
%
% The configuration structure can contain the following
%   cfg.template      = string, filename with the STL model of the sensor sensor or sensor sensor holder
%   cfg.outwardshift  = number, amount to shift the sensor sensors outward from the surface
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.channel       = cell-array, selection of electrode locations at which to place an sensor sensor
%   cfg.rotation      = string, how to determe the cylindrical rotation of the sensor sensor
%   cfg.write         = 'no' or 'yes', write the sensor sensors to STL files
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

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'elec', 'template'});

% set the defaults
cfg.outwardshift  = ft_getopt(cfg, 'outwardshift', 0);
cfg.channel       = ft_getopt(cfg, 'channel', 'all');
cfg.write         = ft_getopt(cfg, 'write', 'no');
cfg.rotation      = ft_getopt(cfg, 'rotation'); % FIXME

% select the channels/electrodes
elec = keepfields(cfg.elec, {'elecpos', 'elecori', 'label'});
cfg.channel = ft_channelselection(cfg.channel, elec.label);
[sel1, sel2] = match_str(cfg.channel, elec.label); % sort them according to cfg.channel
elec.label   = elec.label(sel2);
elec.elecpos = elec.elecpos(sel2,:);
if isfield(elec, 'elecori')
  elec.elecori = elec.elecori(sel2,:);
end
nsens = length(elec.label);

% read the template STL model
template = ft_read_headshape(cfg.template);

% project the electrodes onto the headshape surface
[dum, elec.elecpos] = project_elec(elec.elecpos, headshape.pos, headshape.tri);

if ~isfield(elec, 'elecori')
  % compute the orientation of electrodes
  elec.elecori = normals_elec(elec.elecpos, headshape.pos, headshape.tri);
end % if not elecori

[m, i] = max(headshape.pos(:,3));
reference = headshape.pos(i, :); % approximately the vertex
reference = [0 0 100]; % towards the z-axis

for i=1:nsens
  % shift away from the surface
  t1 = translate([0, 0, cfg.outwardshift]);

  % determine the rotational homogenous matrix
  dirz = elec.elecori(i,:);
  dirz = dirz / norm(dirz);
  dirx = (elec.elecpos(i,:) - reference); % approximate x direction, together with z this defines a plane
  diry = cross(dirz, dirx);
  diry = diry / norm(diry);
  dirx = cross(diry, dirz); % update the actual x direction
  r = [dirx(:) diry(:) dirz(:)];
  r(4,4) = 1;

  % determine the translational homogenous matrix
  t2 = translate(elec.elecpos(i,:));

  % first rotate, then translate
  sensor(i) = ft_transform_geometry(t2 * r * t1, template);

end

if istrue(cfg.write)
  for i=1:nsens
    [p, f, x] = fileparts(cfg.template);
    filename = [f '_' elec.label{i} '.stl'];
    ft_info('writing %s', filename)
    ft_write_headshape(filename, sensor(i), 'format', 'stl');
  end
end

ft_postamble debug
