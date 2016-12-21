function montage = megplanar_orig(cfg, grad)

% This is the original method from Ole.  It has a different way of
% making the coordinate transformation that I do not fully understand.

% Copyright (C) 2004, Ole Jensen
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

neighbsel = cfg.neighbsel;
distance  = cfg.distance;

lab   = grad.label;
tmp   = ft_channelselection('MEG', lab);
sel   = match_str(lab, tmp);
pnt   = grad.chanpos(sel,:);
ori   = grad.chanori(sel,:);
lab   = lab(sel);
Ngrad = length(lab);

gradH = zeros(Ngrad, Ngrad);
gradV = zeros(Ngrad, Ngrad);

% get the positions of bottom and top coil
X  = pnt(:, 1);
Y  = pnt(:, 2);
Z  = pnt(:, 3);
X2 = pnt(:, 1) + ori(:,1);
Y2 = pnt(:, 2) + ori(:,2);
Z2 = pnt(:, 3) + ori(:,3);

for k=1:Ngrad
  % translate the current coil to the origin
  Xc = X - X(k);
  Yc = Y - Y(k);
  Zc = Z - Z(k);
  X2c = X2 - X(k);
  Y2c = Y2 - Y(k);
  Z2c = Z2 - Z(k);

  X = Xc;
  Y = Yc;
  Z = Zc;
  X2 = X2c;
  Y2 = Y2c;
  Z2 = Z2c;

  % rotate around z-axis
  PhiZ = -atan(Y2(k)/(0.00000001+X2(k)));
  Xc  = X*cos(PhiZ) - Y*sin(PhiZ);
  Yc  = X*sin(PhiZ) + Y*cos(PhiZ);
  Zc  = Z;
  X2c = X2*cos(PhiZ) - Y2*sin(PhiZ);
  Y2c = X2*sin(PhiZ) + Y2*cos(PhiZ);
  Z2c = Z2;

  X = Xc;
  Y = Yc;
  Z = Zc;
  X2 = X2c;
  Y2 = Y2c;
  Z2 = Z2c;

  % rotate around y-axis
  PhiY = atan(Z2(k)/(0.00000001+X2(k)));
  Xc  = X*cos(PhiY) + Z*sin(PhiY);
  Yc  = Y;
  Zc  = -X*sin(PhiY) + Z*cos(PhiY);
  X2c = X2*cos(PhiY) + Z2*sin(PhiY);
  Y2c = Y2;
  Z2c = -X2*sin(PhiY) + Z2*cos(PhiY);

  X = Xc;
  Y = Yc;
  Z = Zc;
  X2 = X2c;
  Y2 = Y2c;
  Z2 = Z2c;

  % compute planar gradients in y and z directions
  for l=find(neighbsel(k,:))
    ac = Y(l)/distance(k,l);
    gradH(l,k) = ac/distance(k,l);
    gradH(l,l) = gradH(l,l) - ac/distance(k,l);

    ac = Z(l)/distance(k,l);
    gradV(l,k) = ac/distance(k,l);
    gradV(l,l) = gradV(l,l) - ac/distance(k,l);
  end % for l
end % for k


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
montage.labelorg = lab(:)';                       % describes the columns
montage.tra      = cat(1, gradH, gradV);          % this is the linear projection matrix

  
