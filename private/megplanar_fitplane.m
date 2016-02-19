function montage = megplanar_fitplane(cfg, grad)

% Fit a plane through the B=f(x,y) plane and compute its two gradients
% The first point in the plane is the gradiometer itself,
% the neighbours are the subsequent points. This method also returns the
% offset of the B-plane at each sensor, which is appriximately equal to the
% field itself.

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
% ensure the correct sensor order
[chansel, sel] = match_str(cfg.channel, lab);
pnt   = grad.chanpos(sel,:);
ori   = grad.chanori(sel,:);
lab   = lab(sel);

% we need to ensure that this one is in cfg.channel order - this is done in
% ft_megplanar!
neighbsel = cfg.neighbsel(chansel, chansel);

Ngrad = length(lab);

gradH = zeros(Ngrad, Ngrad);
gradV = zeros(Ngrad, Ngrad);
gradC = zeros(Ngrad, Ngrad);

for chan=1:Ngrad
  % Attach a local coordinate system to this gradiometer:
  % the origin at the location of its bottom coil
  % the z-axis pointing outwards from the head
  % the x-axis pointing horizontal w.r.t. the head
  % the y-axis pointing vertical, i.e. approximately towards the vertex
  this_o = pnt(chan,:);
  this_z = ori(chan,:);
  this_z = this_z / norm(this_z);
  this_x = cross([0 0 1], this_z);
  if all(this_x==0)
    this_x = [1 0 0];
  else
    this_x = this_x / norm(this_x);
  end
  this_y = cross(this_z, this_x);

  % add the relative position in local coordinates to the list of points
  % starting with the channel position itself (which is the local origin)
  x = 0;
  y = 0;
  e = 1;
  neighbindx = find(neighbsel(chan, :));
  for neighb=neighbindx
    vec = pnt(neighb,:) - this_o;          % vector from sensor to neighbour
    x(end+1) = dot(vec, this_x);                % projection along x-axis (horizontal)
    y(end+1) = dot(vec, this_y);                % projection along y-axiz (vertical)
    e(end+1) = 1;                               % required to fit the constant
  end
  A = pinv([x(:) y(:) e(:)]);
  gradH(chan,[chan neighbindx]) = A(1,:);
  gradV(chan,[chan neighbindx]) = A(2,:);
  gradC(chan,[chan neighbindx]) = A(3,:);
end

% rename the labels to match the new channel content
labelH = cell(1, length(lab));
labelV = cell(1, length(lab));
labelC = cell(1, length(lab));
for i=1:length(lab)
  labelH{i} = sprintf('%s_dH', lab{i});
end
for i=1:length(lab)
  labelV{i} = sprintf('%s_dV', lab{i});
end
for i=1:length(lab)
  labelC{i} = sprintf('%s_dC', lab{i});
end

% construct a montage, i.e. a simple linear projection matrix
montage = [];
montage.labelnew = cat(1, labelH(:), labelV(:), labelC(:));   % describes the rows
montage.labelorg = lab(:)';                                   % describes the columns
montage.tra      = cat(1, gradH, gradV, gradC);               % this is the linear projection matrix
