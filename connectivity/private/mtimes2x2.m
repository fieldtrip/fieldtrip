function z = mtimes2x2(x, y)

% MTIMES2X2 compute x*y where the dimensionatity is 2x2xN or 2x2xNxM

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

z     = complex(zeros(size(x)));
xconj = conj(x);

z(1,1,:,:) = x(1,1,:,:).*y(1,1,:,:) + x(1,2,:,:).*y(2,1,:,:);
z(1,2,:,:) = x(1,1,:,:).*y(1,2,:,:) + x(1,2,:,:).*y(2,2,:,:);
z(2,1,:,:) = x(2,1,:,:).*y(1,1,:,:) + x(2,2,:,:).*y(2,1,:,:);
z(2,2,:,:) = x(2,1,:,:).*y(1,2,:,:) + x(2,2,:,:).*y(2,2,:,:);
