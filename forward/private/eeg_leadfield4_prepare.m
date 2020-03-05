function [lut_t, cuf_t] = eeg_leadfield4_prepare(vol, order)

% EEG_LEADFIELD4_PREPARE computes constant factors for series expansion
% for the 4 concentric sphere electric leadfield computation. Calling
% this function speeds up subsequent computations, as the constant
% factors "t" do not have to be computed each time in eeg_leadfield4.
%
% Use as
%   vol.t = eeg_leadfield4_prepare(vol, order);
% where
%   vol.r      radius of the 4 spheres 
%   vol.cond   conductivity of the 4 spheres
% and N is the number of terms for the series (default 60). 
%
% The center of the spheres should be at the origin.
%
% This implementation is adapted from
%   Lutkenhoner, Habilschrift 1992.
% which again is taken from
%   B. N. Cuffin and D. Cohen. Comparion of the Magnetoencephalogram and the Electroencephalogram. Electroencephalogr Clin Neurophysiol, 47:131-146, 1979.
%
% See also EEG_LEADFIELD4

% Copyright (C) 2002, Robert Oostenveld
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

% sort the spheres from the smallest to the largest
[vol.r, indx] = sort(vol.r);
vol.cond      = vol.cond(indx);

r1 = vol.r(1); c1 = vol.cond(1);
r2 = vol.r(2); c2 = vol.cond(2);
r3 = vol.r(3); c3 = vol.cond(3);
r4 = vol.r(4); c4 = vol.cond(4);

if nargin==1 || isempty(order)
  order = 60;
end

% these are the constants of cuffin1979
k1 = c1/c2;
k2 = c2/c3;
k3 = c3/c4;

for n=1:order
  % according to lutkenhoner1992 the constant C is
  % lut_t(n) = ((n*c1/c2+n+1)*(n*c2/c3+n+1)+n*(n+1)*(c1/c2-1)*(c2/c3-1)*(r1/r2)^(2*n+1)) * ...
  %        ((n*c3/c4+n+1)+(n+1)*(c3/c4-1)*(r3/r4)^(2*n+1)) + ...
  %        ((c1/c2-1)*((n+1)*c2/c3+n)*(r1/r3)^(2*n+1)+(n*c1/c2+n+1)*(c2/c3-1)*(r2/r3)^(2*n+1)) * ...
  %        (n+1)*(n*(c3/c4-1)+((n+1)*c3/c4+n)*(r3/r4)^(2*n+1));
  % which can be rewritten as
  lut_t(n) = ((n*k1+n+1)*(n*k2+n+1)+n*(n+1)*(k1-1)*(k2-1)*(r1/r2)^(2*n+1)) * ...
         ((n*k3+n+1)+(n+1)*(k3-1)*(r3/r4)^(2*n+1)) + ...
         ((k1-1)*((n+1)*k2+n)*(r1/r3)^(2*n+1)+(n*k1+n+1)*(k2-1)*(r2/r3)^(2*n+1)) * ...
         (n+1)*(n*(k3-1)+((n+1)*k3+n)*(r3/r4)^(2*n+1));

end

% for debugging purposes, it can also give the constants of cuffin19979
if nargout>1
  % some extra constants of cuffin1979
  b = r1/r4;
  c = r2/r4;
  d = r3/r4;

  % according to cuffin1979 the constant Tau is (re-entered on 25 sept 2002)
  % but this requires also slightly other constants in the eeg_leadfield4 function
  for n=1:order
    cuf_t(n) = d^(2*n+1) * (b^(2*n+1)*n*(k1-1)*(k2-1)*(n+1)...
           + c^(2*n+1)*(k1*n+n+1)*(k2*n+n+1))...
           *((k3*n+n+1)+(n+1)*(k3-1)*d^(2*n+1))...
           +(n+1)*c^(2*n+1)*(b^(2*n+1)*(k1-1)*(k2*n+k2+n)...
           +c^(2*n+1)*(k1*n+n+1)*(k2-1))...
           *(n*(k3-1)+(k3*n+k3+n)*d^(2*n+1));
  end
end
