function [lf] = eeg_infinite_monopole(rd, elc, vol)

% EEG_HALFSPACE_MONOPOLE calculate the halfspace medium leadfield
% on positions pnt for a monopole at position rd and conductivity cond
% The halfspace solution requires a plane dividing a conductive zone of
% conductivity cond, from a non coductive zone (cond = 0)
%       
% [lf] = eeg_halfspace_monopole(rd, elc, cond)
%
% Implemented from Malmivuo J, Plonsey R, Bioelectromagnetism (1993)
% http://www.bem.fi/book/index.htm

% Copyright (C) 2011, Cristiano Micheli
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

siz = size(rd);
if any(siz==1)
  % positions are specified as a single vector
  Npoles = prod(siz)/3;
  rd = rd(:)'; % ensure that it is a row vector
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Npoles = siz(1);
  rd = rd';
  rd = rd(:)'; % ensure that it is a row vector
else
  ft_error('incorrect specification of dipole locations');
end

Nelc     = size(elc,1);
lf       = zeros(Nelc,Npoles); 

for i=1:Npoles
  % this is the position of dipole "i"
  pole1 = rd((1:3) + 3*(i-1));
  
  % distances electrodes - corrent poles
  r1 = elc - ones(Nelc,1) * pole1;
  
  % denominator
  R1 =  (4*pi*1) * (sum(r1' .^2 ) )';
  
  lf(:,i) = (1 ./ R1);

end
