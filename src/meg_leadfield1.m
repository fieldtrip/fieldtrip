function [lf] = meg_leadfield1(dippos, coilpos, coilori)

% MEG_LEADFIELD1 magnetic leadfield for a dipole in a homogenous sphere
%
% Use as
%   [lf] = meg_leadfield1(dippos, coilpos, coilori)
% with input arguments
%   dippos    position of the dipole dipole
%   coilpos   position of the magnetometers
%   coilori   orientation of the magnetometers
%
% The center of the homogenous sphere is in the origin, the field
% of the dipole is not dependent on the sphere radius.
%
% This function is also implemented as MEX file.

% adapted from Luetkenhoener, Habilschrift '92
% optimized for speed using temporary variables
% the mex implementation is a literal copy of this

% Copyright (C) 2002-2008, Robert Oostenveld
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

persistent warning_once
if isempty(warning_once) || ~warning_once
  % the mex file is many times faster than the matlab implementation, hence that is preferred
  % but now we use the matlab implementation as a fallback
  warning_once = true;
  warning('Could not locate the MEX file "%s.%s"', mfilename, mexext);
end

Nchans = size(coilpos, 1);
lf = zeros(Nchans,3);
tmp2 = norm(dippos);

for i=1:Nchans
  r = coilpos(i,:);
  u = coilori(i,:);

  tmp1 = norm(r);
  % tmp2 = norm(R);
  tmp3 = norm(r-dippos);
  tmp4 = dot(r,dippos);
  tmp5 = dot(r,r-dippos);
  tmp6 = dot(dippos,r-dippos);
  tmp7 = (tmp1*tmp2)^2 - tmp4^2;  % cross(r,R)^2

  alpha = 1 / (-tmp3 * (tmp1*tmp3+tmp5));
  A = 1/tmp3 - 2*alpha*tmp2^2 - 1/tmp1;
  B = 2*alpha*tmp4;
  C = -tmp6/(tmp3^3);

  if tmp7<eps
    beta = 0;
  else
    beta = dot(A*r + B*dippos + C*(r-dippos), u)/tmp7;
  end

  lf(i,:) = cross(alpha*u  + beta*r, dippos);
end
lf = 1e-7*lf; % multiply with u0/4pi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product
function [c] = cross(a,b)
c = [a(2)*b(3)-a(3)*b(2) a(3)*b(1)-a(1)*b(3) a(1)*b(2)-a(2)*b(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast dot product
function [c] = dot(a,b)
c = sum(a.*b);

