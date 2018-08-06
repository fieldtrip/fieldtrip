function var = getglobal

% GETGLOBAL gets all global variables and puts them in a structure
%
% Use as
%   var = getglobal;
%   setglobal(var);
%
% See also SETGLOBAL

% Copyright (C) 2011-2012, Robert Oostenveld
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

list = whos('global');

var = [];
for i=1:length(list)
  eval(sprintf('global %s', list(i).name));
  var.(list(i).name) = eval(list(i).name);
end
