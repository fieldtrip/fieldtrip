function x = rmfield(x, key)

% RMFIELD Removes specified field from a CONFIGURATION object.

% Copyright (C) 2012-2015, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if isa(key, 'cell')
  for i=1:numel(key)
    x = rmfield(x, key{i});
  end
else
  x.value     = rmfield(x.value    , key);
  x.assign    = rmfield(x.assign   , key);
  x.reference = rmfield(x.reference, key);
  x.original  = rmfield(x.original , key);
end
