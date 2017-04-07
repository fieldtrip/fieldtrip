function x = set(x, key, val)

% SET Assign a new value to the field of a config object.

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

if ~isfield(x.value, key)
  % initialize the counters for this field
  % see the explaination about side effects of the increment function in config.m
  x.assign.(key)    = deepcopy(0); % ensure that a unique scalar is created for each counter
  x.reference.(key) = deepcopy(0); % ensure that a unique scalar is created for each counter
  x.original.(key)  = deepcopy(0); % ensure that a unique scalar is created for each counter
end

x.value.(key) = val;
increment(x.assign.(key));
