function [elc, lab] = equidistant_locate(pos, tri, front, back, left, right, vertex, numelec, nummidline, numsideline, maxiter, minchange, feedback)

% EQUIDISTANT_LOCATE determines electrode positions that are distributed
% equidistantly on a scalp surface that is described by a triangulation
%
% See also ELEC1020_LOCATE, FT_ELECTRODEPLACEMENT

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

if nargin<12
  feedback = false;
end

if feedback
  figure
end

% the midline is the line connecting the front and the back, over the vertex
% the sideline is the line connecting the front and the back, over the left or right

if isempty(nummidline)
  % FIXME this is not persee a good estimate of the desired number of midline electrodes
  nummidline = floor(sqrt(numelec));
end

if isempty(numsideline)
  % FIXME this is not persee a good estimate of the desired number of sideline electrodes
  numsideline = floor(sqrt(numelec));
end

% to avoid confusion between electrode and headshape positions
headshape.pos = pos;
headshape.tri = tri;
clear pos tri

% determine the surface normals
headshape.nrm = surface_normals(headshape.pos, headshape.tri);

anterior = front - back;
anterior = anterior / norm(anterior);

hasvertex    = ~isempty(vertex);
hasleftright = ~isempty(left) && ~isempty(right);

if ~hasleftright
  center   = (front + back) / 2;
  headsize = norm(front - back);
  % also make an estimate for the left and right points
  if ~hasvertex
    left  = cross([0 0 1], anterior) * headsize/2;
    right = cross(anterior, [0 0 1]) * headsize/2;
  else
    left  = cross(vertex, anterior) * headsize/2;
    right = cross(anterior, vertex) * headsize/2;
  end
else
  center   = (front + back + left + right) / 4;
  headsize = (norm(left-right) + norm(front - back)) / 2; % this is the diameter
end

% determine the homogenous transformation matrix, the FTG coordinate system is defined as:
%   the origin corresponds with pt1
%   the x-axis is along the line from pt1 to pt2
%   the z-axis is orthogonal to the plane spanned by pt1, pt2 and pt3

pt1 = center;
pt2 = center + anterior; % towards the front
pt3 = center + (left - right) / norm(left - right); % towards left

transform = ft_headcoordinates(pt1, pt2, pt3, [], 'ftg');

% transform the coordinates into a "simple" ALS geometry
headshape = ft_transform_geometry(transform, headshape);
front     = ft_transform_geometry(transform, struct('pos', front));
back      = ft_transform_geometry(transform, struct('pos', back));
left      = ft_transform_geometry(transform, struct('pos', left));
right     = ft_transform_geometry(transform, struct('pos', right));
vertex    = ft_transform_geometry(transform, struct('pos', vertex));

% convert these structures back into a simple 1x3 vector
front  = front.pos;
back   = back.pos;
left   = left.pos;
right  = right.pos;
vertex = vertex.pos;

if hasleftright
  % compute the mean of the left and the flipped right electrode position
  leftright = (left + [1 -1 1].*right)/2;
end

% from now on we work in an ALS coordinate system
anterior = [1 0 0];
left     = [0 1 0];
superior = [0 0 1];

if hasvertex
  nummidline = nummidline - 3; % minus vertex, front and back
else
  nummidline = nummidline - 2; % minus front and back
end

if hasleftright
  numsideline = numsideline - 3; % minus left, front and back
else
  numsideline = numsideline - 2; % minus front and back
end

numfixed = 2; % the front and back
numremaining = numelec - numfixed - nummidline - hasvertex - 2*numsideline - 2*hasleftright;
if mod(numremaining, 2)
  % it must be an even number, symmetric over the left and right quadrant
  nummidline = nummidline - 1;
end
numquadrant = (numelec - numfixed - nummidline - hasvertex - 2*numsideline - 2*hasleftright)/2;

% create some electrodes along the midline
midline = zeros(nummidline, 3);
for i=1:nummidline
  angle = i*pi / (nummidline+1);
  midline(i,:) = (cos(angle) * anterior + sin(angle)*superior) * headsize/2;
end

% create some electrodes along the sideline
sideline = zeros(numsideline, 3);
for i=1:numsideline
  angle = i*pi / (numsideline+1);
  sideline(i,:) = (cos(angle) * anterior + sin(angle)*left) * headsize/2;
end

quadrant = uniformquadrant(numquadrant) * headsize/2;
quadrant(:,2) = abs(quadrant(:,2)); % positive x
quadrant(:,3) = abs(quadrant(:,3)); % positive z
for i=1:size(quadrant,1)
  quadrant(i,:) = quadrant(i,:)/norm(quadrant(i,:)) * headsize/2;
end

% the front and back are always fixed
elc = [front; back];
selfixed = [1 2];

if hasvertex
  elc = cat(1, elc, vertex);
  selvertex = size(elc,1);
else
  selvertex = [];
end

elc = cat(1, elc, midline);
selmidline = (size(elc,1)-nummidline+1):size(elc,1);

if hasleftright
  elc = cat(1, elc, leftright);
  selleftright = size(elc,1);
else
  selleftright = [];
end

elc = cat(1, elc, sideline);
selsideline = (size(elc,1)-numsideline+1):size(elc,1);

elc = cat(1, elc, quadrant);
selquadrant = (size(elc,1)-numquadrant+1):size(elc,1);

% project the electrodes onto the mesh surface
[dum, elc] = project_elec(elc, headshape.pos, headshape.tri);

if feedback
  ft_plot_mesh(headshape, 'facecolor', 'skin', 'edgecolor', 'none', 'axes', 'on');
  ft_plot_mesh(elc, 'vertexsize', 20);
  view(-90, 90);
  ft_headlight
end

scale = inf;
change = inf;
iter = 0;

if isempty(minchange)
  minchange = headsize/2000; % approx 0.1 mm for a 200mm head
end

if isempty(maxiter)
  maxiter = 500;
end

while change>minchange && iter<maxiter
  % take the previous positions and copy then over four quadrants
  pos = [
    elc % fixed + vertex + midline + leftright + sideline + quadrant
    elc(selleftright,:) .* repmat([+1 -1 +1], hasleftright, 1)   % mirror leftright
    elc(selsideline,:)  .* repmat([+1 -1 +1], numsideline,  1)   % right sideline
    elc(selquadrant,:)  .* repmat([+1 -1 +1], numquadrant,  1)   % upper right
    elc(selvertex,:)    .* repmat([+1 +1 -1], hasvertex,    1)   % mirror vertex
    elc(selmidline,:)   .* repmat([+1 +1 -1], nummidline,   1)   % lower midline
    elc(selquadrant,:)  .* repmat([+1 +1 -1], numquadrant,  1)   % lower left
    elc(selquadrant,:)  .* repmat([+1 -1 -1], numquadrant,  1)   % lower right
    ];

  % compute the repulsion force between all points
  shift = zeros(size(elc));
  for i=1:size(elc,1)
    for j=1:size(pos,1)
      if i==j
        continue
      else
        v = pos(i,:) - pos(j,:);
        d = norm(v);
        v = v/d;
        w = 1/d^2;
        shift(i,:) = shift(i,:) + w * v;
      end
    end
  end

  % none of the points is to shift with a too large amount
  mx = 2*median(abs(shift(:,1)));
  my = 2*median(abs(shift(:,2)));
  mz = 2*median(abs(shift(:,3)));
  sel = shift(:,1)>mx; shift(sel,1) = mx;
  sel = shift(:,2)>my; shift(sel,2) = my;
  sel = shift(:,3)>mz; shift(sel,3) = mz;
  sel = shift(:,1)<-mx; shift(sel,1) = -mx;
  sel = shift(:,2)<-my; shift(sel,2) = -my;
  sel = shift(:,3)<-mz; shift(sel,3) = -mz;

  % the scaling is only computed once
  if isinf(scale)
    scale = 3/(mx+my+mz);
  end
  shift = shift * scale;

  % the fixed points are not to shift
  shift(selfixed,:)     = 0;
  shift(selvertex,:)    = 0;
  shift(selleftright,:) = 0;

  % the midline points are not to shift in the y-direction
  shift(selmidline,2) = 0;

  % the sideline points are not to shift in the z-direction
  shift(selsideline,3) = 0;

  % none of the points is to shift below z=0
  sel = (elc(:,3)+shift(:,3))<0;
  shift(sel,3) = 0;

  % update the left quadrant, drop the other three quadrants
  pos = pos(1:size(elc,1),:) + shift;

  % project the updated points onto the headsurface
  [dum, pos] = project_elec(pos, headshape.pos, headshape.tri);

  % these have the tendency to drift away due to repeated projections, keep them fixed
  pos(selmidline,2)   = 0;
  pos(selsideline,3)  = 0;

  % compute the average shift of the points
  change = sum(sqrt(sum((elc - pos).^2,2)));
  iter = iter + 1;

  if feedback && mod(iter,10)==0
    cla
    ft_info off
    ft_plot_mesh(pos, 'axes', 'on', 'vertexsize', 20);
    ft_info on
    drawnow
    ft_info('iter = %d, change = %f\n', iter, change)
  end

  % repeat with the updated positions
  elc = pos;
end % while

% copy the electrode positions from the left to the right quadrant
elc = [
  elc % fixed + vertex + midline + leftright + sideline + quadrant
  elc(selleftright,:) .* repmat([+1 -1 +1], hasleftright, 1)  % mirror leftright
  elc(selsideline,:)  .* repmat([+1 -1 +1], numsideline,  1)  % right sideline
  elc(selquadrant,:)  .* repmat([+1 -1 +1], numquadrant,  1)  % upper right
  ];

% sanity check
assert(size(elc,1)==numelec);

% FIXME this is currently hard-coded
sortorder = 'lrfb';

switch sortorder
  case 'azel'
    % sort them according to the azimuth and elevation
    [az,el] = cart2sph(elc(:,1), elc(:,2), elc(:,3));
    [dum, indx] = sortrows([az el], [2 1], 'descend');
    elc = elc(indx,:);
  case 'lrfb'
    % sort from left-to-right, then from front-to-back
    lr = elc(:,2);
    fb = elc(:,1);
    [dum, indx] = sortrows([lr fb], [2 1], 'descend');
    elc = elc(indx,:);
end

% construct electrode labels
lab = cell(numelec, 1);
for i=1:numelec
  lab{i} = sprintf('%d', i);
end

% transform the electrodes back to the original coordinate system
elc = ft_transform_geometry(inv(transform), struct('pos', elc));

% convert the structures back into a simple Nx3 matrix
elc = elc.pos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = uniformquadrant(n)
% this function distributes N points more or less uniform
% over the left quadrant in 3D space

az = 0:(n-1);
el = 0:(n-1);

m  = floor(sqrt(n));
az = mod(az, m)/(m-1);
el = el / (n-1);

az = pi/(1*m) + az * pi/1 * (1 - 2/m);
el = pi/(2*m) + el * pi/2 * (1 - 2/m);

[x, y, z] = sph2cart(az, el, ones(1,n));

pos = [x' y' z'];