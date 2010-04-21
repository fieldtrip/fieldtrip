function grad = itab2grad(header_info)

% ITAB2GRAD converts the original Chieti ITAB header structure into a gradiometer
% definition that is compatible with FieldTrip forward and inverse computations
%
% See also READ_HEADER

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
for i=1:header_info.nmagch
  grad.label{i}   = header_info.ch(i).label;
  grad.pnt(i,1:3) = [header_info.ch(i).pos(1).r_s.comp];
  grad.ori(i,1:3) = [header_info.ch(i).pos(1).u_s.comp];
end
grad.unit  = 'mm';
grad.tra   = eye(header_info.nmagch);
grad.label = grad.label(:); % should be column vector
