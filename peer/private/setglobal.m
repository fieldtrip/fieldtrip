function setglobal(var)

% SETGLOBAL assigns global variables from a structure
%
% Use as
%   var = getglobal;
%   setglobal(var);
%
% See also GETGLOBAL

% Copyright (C) 2012, Robert Oostenveld
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

if isempty(var)
  return
end

list = fieldnames(var);

for i=1:length(list)
  eval(sprintf('global %s', list{i}));
  eval(sprintf('%s = var.%s;', list{i}, list{i}));
end
