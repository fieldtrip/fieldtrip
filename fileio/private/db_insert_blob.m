function db_insert_blob(tablename, fieldname, s)

% DB_INSERT_BLOB converts a Matlab variable of arbitrary type into
% a binary stream and inserts this stream into a binary blob in the
% database table.
%
% Use as
%   db_insert_blob(tablename, fieldname, s)
%
% See also DB_OPEN, DB_SELECT, DB_SELECT_BLOB, DB_INSERT, DB_CLOSE

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
cmd    = sprintf('insert into %s ( %s ) values ( "{M}" )', tablename, fieldname);
blob   = struct2msg(s, format);

% execute the query, this requires the mym version of the mysql mex file
mysql(str, blob);
