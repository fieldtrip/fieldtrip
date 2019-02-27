function [lf] = eeg_infinite_monopole(monpos, elc, vol)

% EEG_INFINITE_MONOPOLE calculate the infinite medium potential for a monopole
%
% Use as
%   [lf] = eeg_infinite_monopole(monpos, elc, vol)
%
% Implemented from Malmivuo J, Plonsey R, Bioelectromagnetism (1993)
% http://www.bem.fi/book/08/08.htm
%
% See also EEG_INFINITE_DIPOLE, EEG_HALFSPACE_DIPOLE, EEG_HALFSPACE_MONOPOLE

% Copyright (C) 2011, Cristiano Micheli
% Copyright (C) 2019, Robert Oostenveld
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

if ~isstruct(vol)
  % it only represents the conductivity, make a structure out of it
  vol = struct('cond', vol);
end

siz = size(monpos);
if any(siz==1)
  % positions are specified as a single vector
  Npoles = prod(siz)/3;
  monpos = monpos(:)'; % ensure that it is a row vector
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Npoles = siz(1);
  monpos = monpos';
  monpos = monpos(:)'; % ensure that it is a row vector
else
  ft_error('incorrect specification of monopole locations');
end

cond     = vol.cond;
Nelc     = size(elc,1);
lf       = zeros(Nelc,Npoles);

mu0   = 4*pi*1e-7;         % Permeability of free space
c     = 2.99792458 * 1e8;  % Speed of light
e0    = 1 / (mu0*c^2);     % Permittivity of Free Space

for i=1:Npoles
  % this is the position of monopole "i"
  monopole = monpos((1:3) + 3*(i-1));
  
  % distances from electrodes to monopole
  r = elc - ones(Nelc,1) * monopole;
  r = sqrt(sum(r.^2,2));
  
  lf(:,i) = 1 ./ (4*pi*cond*r);
end
