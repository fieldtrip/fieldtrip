function A = reorderdim(A, dim, inds)

% REORDERDIM reorders array A along dimension dim with the specified
% indices inds. The following should output 1:
%
%   B1 = reorderdim(A,2,[1 3 2]);
%   B2 = A(:,[1 3 2],:,:);
%
%   all(B1(:) == B2(:))
%
% The main use for this function is when a selection as displayed above
% needs to be made when the number of dimensions of A is only known at
% runtime and not at 'code'-time (i.e. when A can have arbitrary
% dimensions).
%

% Copyright (C) 2013 Eelke Spaak
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

% the functionality is achieved through some ugly string manipulation

indString = 'A=A(';
for l = 1:numel(size(A))
  if l == dim
    indString = [indString mat2str(inds), ','];
  else
    indString = [indString ':,'];
  end
end
indString = [indString(1:end-1) ');'];
eval(indString);

end
