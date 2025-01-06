function s = db_select_blob(tablename, fieldname, sel)

% DB_SELECT_BLOB selects a binary blob from a database table and converts
% it back into a Matlab variable. The variable can be of an arbitrary type.
%
% Use as
%   s = db_select_blob(tablename, fieldname)
%   s = db_select_blob(tablename, fieldname, num)
%
% The optional argument num allows you to select a specific row number.
%
% See also DB_OPEN, DB_INSERT, DB_SELECT, DB_INSERT_BLOB, DB_CLOSE

% Copyright (C) 2007, Robert Oostenveld
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

format = 'serialize';

if nargin>2
  cmd = sprintf('SELECT %s FROM %s LIMIT %d,1', fieldname, tablename, sel-1);
else
  cmd = sprintf('SELECT %s FROM %s', fieldname, tablename);
end

% execute the query, this requires the mym version of the mysql mex file
s = mysql(cmd);

if numel(s)==1
  s = msg2struct(s{1}, format);
else
  for i=1:numel(s)
    s{i} = msg2struct(s{i}, format);
  end
end
