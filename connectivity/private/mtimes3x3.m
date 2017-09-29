function z = mtimes3x3(x, y)

% MTIMES3X3 compute x*y where the dimensionatity is 3x3xN or 3x3xNxM

% Copyright (C) 2017, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

z     = complex(zeros(size(x)));
%xconj = conj(x);

z(1,1,:,:) = x(1,1,:,:).*y(1,1,:,:) + x(1,2,:,:).*y(2,1,:,:) + x(1,3,:,:).*y(3,1,:,:);
z(1,2,:,:) = x(1,1,:,:).*y(1,2,:,:) + x(1,2,:,:).*y(2,2,:,:) + x(1,3,:,:).*y(3,2,:,:);
z(1,3,:,:) = x(1,1,:,:).*y(1,3,:,:) + x(1,2,:,:).*y(2,3,:,:) + x(1,3,:,:).*y(3,3,:,:);
z(2,1,:,:) = x(2,1,:,:).*y(1,1,:,:) + x(2,2,:,:).*y(2,1,:,:) + x(2,3,:,:).*y(3,1,:,:);
z(2,2,:,:) = x(2,1,:,:).*y(1,2,:,:) + x(2,2,:,:).*y(2,2,:,:) + x(2,3,:,:).*y(3,2,:,:);
z(2,3,:,:) = x(2,1,:,:).*y(1,3,:,:) + x(2,2,:,:).*y(2,3,:,:) + x(2,3,:,:).*y(3,3,:,:);
z(3,1,:,:) = x(3,1,:,:).*y(1,1,:,:) + x(3,2,:,:).*y(2,1,:,:) + x(3,3,:,:).*y(3,1,:,:);
z(3,2,:,:) = x(3,1,:,:).*y(1,2,:,:) + x(3,2,:,:).*y(2,2,:,:) + x(3,3,:,:).*y(3,2,:,:);
z(3,3,:,:) = x(3,1,:,:).*y(1,3,:,:) + x(3,2,:,:).*y(2,3,:,:) + x(3,3,:,:).*y(3,3,:,:);
