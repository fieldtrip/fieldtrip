function z = sandwich2x2(x, y)

% SANDWICH2X2 compute x*y*x' provided y is Hermitian and dimensionality is 2x2xN

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

% FIXME build in check for hermitianity
z     = complex(zeros(size(x)));
xconj = conj(x);
xabs2 = abs(x).^2;

z(1,1,:,:) = xabs2(1,1,:,:) .* y(1,1,:,:) + ...
             xabs2(1,2,:,:) .* y(2,2,:,:) + ...
           2.*real(y(2,1,:,:).*xconj(1,1,:,:).*x(1,2,:,:));
z(2,1,:,:) = y(1,1,:,:).*xconj(1,1,:,:).*x(2,1,:,:) + ...
           y(2,1,:,:).*xconj(1,1,:,:).*x(2,2,:,:) + ...
           y(1,2,:,:).*xconj(1,2,:,:).*x(2,1,:,:) + ...
           y(2,2,:,:).*xconj(1,2,:,:).*x(2,2,:,:);
z(1,2,:,:) = conj(z(2,1,:,:));
z(2,2,:,:) = xabs2(2,1,:,:) .* y(1,1,:,:) + ...
             xabs2(2,2,:,:) .* y(2,2,:,:) + ...
           2.*real(y(1,2,:,:).*xconj(2,2,:,:).*x(2,1,:,:));

%b1 b2     a1 a2'   b1' b3'
%b3 b4     a2 a3    b2' b4'
%
%b1*a1+b2*a2  b1*a2'+b2*a3  b1' b3'
%b3*a1+b4*a2  b3*a2'+b4*a3  b2' b4'
%
%b1*a1*b1'+b2*a2*b1'+b1*a2'*b2'+b2*a3*b2' b1*a1*b3'+b2*a2*b3'+b1*a2'*b4'+b2*a3*b4'
%b3*a1*b1'+b4*a2*b1'+b3*a2'*b2'+b4*a3*b2' b3*a1*b3'+b4*a2*b3'+b3*a2'*b4'+b4*a3*b4'
%
%a1*abs(b1)^2 + a2*(b1'*b2) + a2'*(b1*b2') + a3*abs(b2)^2    a1*b1*b3'    + a2*b2*b3'   + a2'*b1*b4'   + a3*b2*b4'
%a1*b1'*b3    + a2*b1'*b4   + a2'*b2'*b3   + a3*b2'*b4       a1*abs(b3)^2 + a2*(b3'*b4) + a2'*(b3*b4') + a3*abs(b4)^2
