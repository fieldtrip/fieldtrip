function z = sandwich3x3(x, y)

% SANDWICH3X3 compute x*y*x' provided y is Hermitian and dimensionality is 3x3xN

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

% FIXME build in check for hermitianity
z     = complex(zeros(size(x)));
xconj = conj(x);
xabs2 = abs(x).^2;

z(1,1,:,:) = xabs2(1,1,:,:).*y(1,1,:,:)  + ...
             xabs2(1,2,:,:).*y(2,2,:,:)  + ...
             xabs2(1,3,:,:).*y(3,3,:,:)  + ...
             2.*real( x(1,2,:,:).*y(2,1,:,:).*xconj(1,1,:,:) ) + ...
             2.*real( x(1,3,:,:).*y(3,1,:,:).*xconj(1,1,:,:) ) + ...
             2.*real( x(1,3,:,:).*y(3,2,:,:).*xconj(1,2,:,:) );
             
             %(x(1,2,:,:).*y(2,1,:,:) + x(1,3,:,:).*y(3,1,:,:)).*xconj(1,1,:,:) + ...
             %(x(1,1,:,:).*y(1,2,:,:) + x(1,3,:,:).*y(3,2,:,:)).*xconj(1,2,:,:) + ...
             %(x(1,1,:,:).*y(1,3,:,:) + x(1,2,:,:).*y(2,3,:,:)).*xconj(1,3,:,:);

z(2,1,:,:) = (x(2,1,:,:).*y(1,1,:,:) + x(2,2,:,:).*y(2,1,:,:) + x(2,3,:,:).*y(3,1,:,:)).*xconj(1,1,:,:) + ...
             (x(2,1,:,:).*y(1,2,:,:) + x(2,2,:,:).*y(2,2,:,:) + x(2,3,:,:).*y(3,2,:,:)).*xconj(1,2,:,:) + ...
             (x(2,1,:,:).*y(1,3,:,:) + x(2,2,:,:).*y(2,3,:,:) + x(2,3,:,:).*y(3,3,:,:)).*xconj(1,3,:,:);

z(3,1,:,:) = (x(3,1,:,:).*y(1,1,:,:) + x(3,2,:,:).*y(2,1,:,:) + x(3,3,:,:).*y(3,1,:,:)).*xconj(1,1,:,:) + ...
             (x(3,1,:,:).*y(1,2,:,:) + x(3,2,:,:).*y(2,2,:,:) + x(3,3,:,:).*y(3,2,:,:)).*xconj(1,2,:,:) + ...
             (x(3,1,:,:).*y(1,3,:,:) + x(3,2,:,:).*y(2,3,:,:) + x(3,3,:,:).*y(3,3,:,:)).*xconj(1,3,:,:);

z(1,2,:,:) = conj(z(2,1,:,:));

z(2,2,:,:) = xabs2(2,1,:,:).*y(1,1,:,:)  + ...
             xabs2(2,2,:,:).*y(2,2,:,:)  + ...
             xabs2(2,3,:,:).*y(3,3,:,:)  + ...
             2.*real( x(2,2,:,:).*y(2,1,:,:).*xconj(2,1,:,:) ) + ...
             2.*real( x(2,3,:,:).*y(3,1,:,:).*xconj(2,1,:,:) ) + ...
             2.*real( x(2,3,:,:).*y(3,2,:,:).*xconj(2,2,:,:) );
                           
             %(x(2,2,:,:).*y(2,1,:,:) + x(2,3,:,:).*y(3,1,:,:)).*xconj(2,1,:,:) + ...
             %(x(2,1,:,:).*y(1,2,:,:) + x(2,3,:,:).*y(3,2,:,:)).*xconj(2,2,:,:) + ...
             %(x(2,1,:,:).*y(1,3,:,:) + x(2,2,:,:).*y(2,3,:,:)).*xconj(2,3,:,:);

z(3,2,:,:) = (x(3,1,:,:).*y(1,1,:,:) + x(3,2,:,:).*y(2,1,:,:) + x(3,3,:,:).*y(3,1,:,:)).*xconj(2,1,:,:) + ...
             (x(3,1,:,:).*y(1,2,:,:) + x(3,2,:,:).*y(2,2,:,:) + x(3,3,:,:).*y(3,2,:,:)).*xconj(2,2,:,:) + ...
             (x(3,1,:,:).*y(1,3,:,:) + x(3,2,:,:).*y(2,3,:,:) + x(3,3,:,:).*y(3,3,:,:)).*xconj(2,3,:,:);

z(1,3,:,:) = conj(z(3,1,:,:));

z(2,3,:,:) = conj(z(3,2,:,:));

z(3,3,:,:) = xabs2(3,1,:,:).*y(1,1,:,:)  + ...
             xabs2(3,2,:,:).*y(2,2,:,:)  + ...
             xabs2(3,3,:,:).*y(3,3,:,:)  + ...
             2.*real( x(3,2,:,:).*y(2,1,:,:).*xconj(3,1,:,:) ) + ...
             2.*real( x(3,3,:,:).*y(3,1,:,:).*xconj(3,1,:,:) ) + ...
             2.*real( x(3,3,:,:).*y(3,2,:,:).*xconj(3,2,:,:) );
            
             %(x(3,2,:,:).*y(2,1,:,:) + x(3,3,:,:).*y(3,1,:,:)).*xconj(3,1,:,:) + ...
             %(x(3,1,:,:).*y(1,2,:,:) + x(3,3,:,:).*y(3,2,:,:)).*xconj(3,2,:,:) + ...
             %(x(3,1,:,:).*y(1,3,:,:) + x(3,2,:,:).*y(2,3,:,:)).*xconj(3,3,:,:);
%
%b1 b2 b7    a1 a2' a4'   b1' b3' b5'
%b3 b4 b8    a2 a3  a5'   b2' b4' b6'
%b5 b6 b9    a4 a5  a6    b7' b8' b9'
%
%b1*a1+b2*a2+b7*a4  b1*a2'+b2*a3+b7*a5  b1*a4'+b2*a5'+b7*a6   b1' b3' b5'
%b3*a1+b4*a2+b8*a4  b3*a2'+b4*a3+b8*a5  b3*a4'+b4*a5'+b8*a6   b2' b4' b6'
%b5*a1+b6*a2+b9*a4  b5*a2'+b6*a3+b9*a5  b5*a4'+b6*a5'+b9*a6   b7' b8' b9'
%
%

%b1*a1*b1'+b2*a2*b1'+b1*a2'*b2'+b2*a3*b2' b1*a1*b3'+b2*a2*b3'+b1*a2'*b4'+b2*a3*b4'
%b3*a1*b1'+b4*a2*b1'+b3*a2'*b2'+b4*a3*b2' b3*a1*b3'+b4*a2*b3'+b3*a2'*b4'+b4*a3*b4'
%
%a1*abs(b1)^2 + a2*(b1'*b2) + a2'*(b1*b2') + a3*abs(b2)^2    a1*b1*b3'    + a2*b2*b3'   + a2'*b1*b4'   + a3*b2*b4'
%a1*b1'*b3    + a2*b1'*b4   + a2'*b2'*b3   + a3*b2'*b4       a1*abs(b3)^2 + a2*(b3'*b4) + a2'*(b3*b4') + a3*abs(b4)^2
