function d = det2x2(x)

% DET2X2 computes determinant of matrix x, using explicit analytic definition
% if size(x,1) < 4, otherwise use MATLAB det-function

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

siz = size(x);
if all(siz(1:2)==2),
  d = x(1,1,:,:).*x(2,2,:,:) - x(1,2,:,:).*x(2,1,:,:);
elseif all(siz(1:2)==3),
  d = x(1,1,:,:).*x(2,2,:,:).*x(3,3,:,:) - ...
      x(1,1,:,:).*x(2,3,:,:).*x(3,2,:,:) - ...
      x(1,2,:,:).*x(2,1,:,:).*x(3,3,:,:) + ...
      x(1,2,:,:).*x(2,3,:,:).*x(3,1,:,:) + ...
      x(1,3,:,:).*x(2,1,:,:).*x(3,2,:,:) - ...
      x(1,3,:,:).*x(2,2,:,:).*x(3,1,:,:);
elseif numel(siz)==2,
  d = det(x);
else
  error('not implemented');
  % write for loop for the higher dimensions, using normal inv
end
