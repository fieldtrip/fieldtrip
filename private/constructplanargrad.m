function [planar] = constructplanargrad(cfg, grad)

% CONSTRUCTPLANARGRAD constructs a planar gradiometer array from an axial gradiometer
% definition. This can be used to compute the planar field gradient for a known
% (estimated) source configuration.
% 
% Use as
%   [grad_planar] = constructplanargrad(cfg, grad_axial)
%
% Where cfg contains the following configuration details
%   cfg.baseline_axial   = number (default is 5)
%   cfg.baseline_planar  = number (default is 0.5)
%   cfg.planaraxial      = 'no' or 'yes' (default)
% 
% The option planaraxial='yes' specifies that the planar gradiometers
% should consist of axial gradiometers, to make them comparable with
% Ole Jensens planar gradient computation. If planaraxial='no', the
% planar gradiometers will be more or less similar to the Neuromag
% system.
%
% The input grad can be a CTF type axial gradiometer definition, but
% just as well be a magnetometer definition. This function only assumes
% that
%   grad.coilpos
%   grad.coilori
%   grad.label
% exist and that the first Nlabel channels in pnt and ori should be
% used to compute the position of the coils in the planar gradiometer
% channels.

% Copyright (C) 2004, Robert Oostenveld
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

if ~isfield(cfg, 'planaraxial'),     cfg.planaraxial = 'yes';   end
if ~isfield(cfg, 'baseline_axial'),  cfg.baseline_axial  = 5;   end
if ~isfield(cfg, 'baseline_planar'), cfg.baseline_planar = 0.5; end

Nchan = length(grad.label);

% these will hold all the coil positions
lo_posx = zeros(Nchan,3);
lo_negx = zeros(Nchan,3);
lo_posy = zeros(Nchan,3);
lo_negy = zeros(Nchan,3);
hi_posx = zeros(Nchan,3);
hi_negx = zeros(Nchan,3);
hi_posy = zeros(Nchan,3);
hi_negy = zeros(Nchan,3);

for chan=1:Nchan
  % Attach a local coordinate system to this gradiometer:
  %   the origin at the location of its bottom coil
  %   the z-axis pointing outwards from the head
  %   the x-axis pointing horizontal w.r.t. the head
  %   the y-axis pointing vertical, i.e. approximately towards the vertex
  this_o = grad.chanpos(chan,:);
  this_z = grad.chanori(chan,:);          
  this_z = this_z / norm(this_z);
  this_x = cross([0 0 1], this_z);
  if all(this_x==0)
    this_x = [1 0 0];
  else
    this_x = this_x / norm(this_x);
  end
  this_y = cross(this_z, this_x);

  % compute the position of all the 8 coils per channel
  lo_posx(chan,:) = this_o + (cfg.baseline_planar/2) * this_x;
  lo_negx(chan,:) = this_o - (cfg.baseline_planar/2) * this_x;
  lo_posy(chan,:) = this_o + (cfg.baseline_planar/2) * this_y;
  lo_negy(chan,:) = this_o - (cfg.baseline_planar/2) * this_y;
  hi_posx(chan,:) = lo_posx(chan,:) + cfg.baseline_axial * this_z;
  hi_negx(chan,:) = lo_negx(chan,:) + cfg.baseline_axial * this_z;
  hi_posy(chan,:) = lo_posy(chan,:) + cfg.baseline_axial * this_z;
  hi_negy(chan,:) = lo_negy(chan,:) + cfg.baseline_axial * this_z;
end

% start with an empty planar gradiometer definition
planar = [];

if strcmp(cfg.planaraxial, 'yes')
  % combine all the 8 coils into a single sensor
  planar.coilpos = [
    lo_posx
    lo_negx
    lo_posy
    lo_negy
    hi_posx
    hi_negx
    hi_posy
    hi_negy
  ];

  % the orientation of all the coils of a single sensor should be the same
  planar.coilori = [
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
  ];

  e = eye(Nchan);
  z = zeros(Nchan);

  % the linear combination matrix should be 2*Nchan x 8*Nchan
  planar.tra = [
    e -e  z  z -e  e  z  z    % this is for the horizontal gradients
    z  z  e -e  z  z -e  e    % this is for the vertical gradients
  ];

else
  % combine only the 4 lower coils into a single sensor
  planar.coilpos = [
    posx
    negx
    posy
    negy
  ];

  % the orientation of all the coils of a single gradiometer should be the same
  planar.coilori = [
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
    grad.coilori(1:Nchan,:)
  ];

  e = eye(Nchan);
  z = zeros(Nchan);

  % the linear combination matrix should be 2*Nchan x 4*Nchan
  planar.tra = [
    e -e  z  z    % this is for the horizontal gradients
    z  z  e -e    % this is for the vertical gradients
  ];

end

for chan=1:Nchan
  planar.label{chan      } = [grad.label{chan} '_dH'];
  planar.label{chan+Nchan} = [grad.label{chan} '_dV'];
end

planar.label = planar.label(:);
planar.tra   = planar.tra / cfg.baseline_planar;
planar.chanpos = [grad.chanpos; grad.chanpos];
planar.chanori = [grad.chanori; grad.chanori];

try
  planar.unit  = grad.unit;
end

% add information about the version of this function to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id   = '$Id$';

% rememember the exact configuration details in the output
planar.cfg = cfg;

