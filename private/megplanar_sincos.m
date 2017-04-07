function montage = megplanar_sincos(cfg, grad)

% This attempts to re-implements Ole's method, exept that the definition of the
% horizontal and vertical direction is different.

% Copyright (C) 2004-2009, Robert Oostenveld
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

lab   = grad.label;
type  = grad.chantype;

% ensure correct order
% cfg.channel       = ft_channelselection(cfg.channel, lab);
[chansel, labsel] = match_str(cfg.channel, lab);
lab               = lab(labsel);
type              = type(labsel);

% we need to ensure that this one is in cfg.channel order - this is done in
% ft_megplanar!
neighbsel         = cfg.neighbsel(chansel, chansel);

% sel   = match_str(lab, tmp);
pnt   = grad.chanpos(labsel,:);
ori   = grad.chanori(labsel,:);
Ngrad = numel(labsel);

gradH = zeros(Ngrad, Ngrad);
gradV = zeros(Ngrad, Ngrad);

for chan=1:Ngrad
  % Attach a local coordinate system to this gradiometer:
  % the origin at the location of its bottom coil
  % the z-axis pointing outwards from the head
  % the x-axis pointing horizontal w.r.t. the head
  % the y-axis pointing vertical, i.e. approximately towards the vertex
  this_o = pnt(chan,:);
  this_z = ori(chan,:);          this_z = this_z / norm(this_z);
  this_x = cross([0 0 1], this_z);
  if all(this_x==0)
    this_x = [1 0 0];
  else
    this_x = this_x / norm(this_x);
  end
  this_y = cross(this_z, this_x);

  for neighb=find(neighbsel(chan, :))
    vec = pnt(neighb,:) - this_o;    % vector from sensor to neighbour
    proj_x = dot(vec, this_x);            % projection along x-axis (horizontal)
    proj_y = dot(vec, this_y);            % projection along y-axiz (vertical)
    proj_z = dot(vec, this_z);            % projection along z-axis

    gradH(chan, chan)   = gradH(chan,chan)    - proj_x / (norm(vec).^2);
    gradH(chan, neighb) =                       proj_x / (norm(vec).^2);
    gradV(chan, chan)   = gradV(chan,chan)    - proj_y / (norm(vec).^2);
    gradV(chan, neighb) =                       proj_y / (norm(vec).^2);
  end
end

% rename the labels to match the new channel content
labelH = cell(1, length(lab));
labelV = cell(1, length(lab));
for i=1:length(lab)
  labelH{i} = sprintf('%s_dH', lab{i});
end
for i=1:length(lab)
  labelV{i} = sprintf('%s_dV', lab{i});
end

% construct a montage, i.e. a simple linear projection matrix
montage = [];
montage.labelnew = cat(1, labelH(:), labelV(:));  % describes the rows
montage.labelold = lab(:)';                       % describes the columns
montage.chantypenew = repmat(type, [2 1]);
montage.chantypeold = type;
montage.tra      = cat(1, gradH, gradV);          % this is the linear projection matrix
