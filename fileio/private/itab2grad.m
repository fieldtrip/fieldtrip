function grad = itab2grad(header_info)

% ITAB2GRAD converts the original Chieti ITAB header structure into
% a gradiometer definition that is compatible with FieldTrip forward
% and inverse computations
%
% See also CTF2GRAD, BTI2GRAD, FIF2GRAD, MNE2GRAD, YOKOGAWA2GRAD,
% FT_READ_SENS, FT_READ_HEADER

% Copyright (C) 2009, Robert Oostenveld, Donders Institute for Brain, Cognition and Behaviour
% Copyright (C) 2009, Stefania Della Penna, ITAB, University Chiety, Italy
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

grad = struct;

chan = 0;
coil = 0;
for i=1:640
  if header_info.ch(i).type==2
    % it is a MAG channel
    chan = chan+1;
    grad.label{chan} = header_info.ch(i).label; 
    for j=1:header_info.ch(i).ncoils
      coil = coil+1;
      grad.coilpos(coil,:) = header_info.ch(i).position(j).r_s;
      grad.coilori(coil,:) = header_info.ch(i).position(j).u_s;
      grad.tra(chan,coil) = header_info.ch(i).wgt(j);
    end    
  else
    % skip all other channels
  end
end

grad.unit  = 'mm';
grad.label = grad.label(:); % should be column vector
