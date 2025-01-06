function db_insert(tablename, s)

% DB_INSERT inserts a structure into a database table. Each field of
% the structure should correspond with one of the fields in the table.
%
% Use as
%   db_insert(tablename, s)
%
% See also DB_OPEN, DB_SELECT, DB_SELECT_BLOB, DB_INSERT_BLOB, DB_CLOSE

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

if isempty(s)
  % the input consists of an empty structure. so nothing to do
  return
end

f  = fieldnames(s);
v = {};
t = {};
l = length(f);
for i=1:l
  v{i} = getfield(s, f{i});
  t{i} = class(v{i});
end

for j=1:numel(s)
  % loop over all elements of a structure array

  cmd = ['INSERT INTO ' tablename ' ( '];

  % these are the fieldnames
  for i=1:l
    cmd = [cmd f{i} ', '];
  end
  cmd = cmd(1:(end-2));  % remove the last ', '
  cmd = [cmd ' ) VALUES ( '];

  % these are the corresponding values
  for i=1:l
    switch t{i}
      case 'char'
        if isempty(v{i})
          cmd = [cmd '""' ];
        else
          cmd = [cmd '"' v{i} '"' ];
        end

      case {'int8', 'int16', 'int32', 'uint8', 'uint16', 'uint32', 'double'}
        if isempty(v{i})
          cmd = [cmd '""' ];
        else
          cmd = [cmd num2str(v{i})];
        end

      otherwise
        ft_error('unsuported data type in structure');
    end
    cmd = [cmd ', '];
  end

  cmd = cmd(1:(end-2));  % remove the last ', '
  cmd = [cmd ' ) ;'];

  % execute the query
  mysql(cmd);

end %for numel(s)

