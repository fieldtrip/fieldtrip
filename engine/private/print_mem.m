function str = print_mem(val)

% PRINT_MEM helper function for pretty-printing memory in GB, MB, ...

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

if val<1024
  str = sprintf('%d bytes', val);
elseif val<1024^2
  str = sprintf('%.1f KB', val/1024);
elseif val<1024^3
  str = sprintf('%.1f MB', val/1024^2);
else
  str = sprintf('%.1f GB', val/1024^3);
end

