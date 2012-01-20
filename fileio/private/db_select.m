function s = db_select(tablename, fields, sel)

% DB_SELECT selects data from a database table and converts it into a
% Matlab structure.Each of the fields in the database table will be
% represented as field in the strucure.
%
% Use as
%   s = db_select(tablename, fields)
%   s = db_select(tablename, fields, num)
%
% The optional argument num allows you to select a specific row number.

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: db_select.m,v $
% Revision 1.1  2007/11/07 10:50:46  roboos
% created helper functions for easy access to a MySQL database table using a structure for representing the data
%

nfields  = length(fields);
fieldstr = sprintf('%s, ', fields{:});
fieldstr = fieldstr(1:(end-2)); % remove the last comma

if nargin>2
  cmd = sprintf('SELECT %s FROM %s LIMIT %d,1', fieldstr, tablename, sel-1);
else
  cmd = sprintf('SELECT %s FROM %s', fieldstr, tablename);
end

% execute the query
switch nfields
  case 1
    [t{1}] = mysql(cmd);
  case 2
    [t{1}, t{2}] = mysql(cmd);
  case 3
    [t{1}, t{2}, t{3}] = mysql(cmd);
  case 4
    [t{1}, t{2}, t{3}, t{4}] = mysql(cmd);
  case 5
    [t{1}, t{2}, t{3}, t{4}, t{5}] = mysql(cmd);
  case 6
    [t{1}, t{2}, t{3}, t{4}, t{5}, t{6}] = mysql(cmd);
  case 7
    [t{1}, t{2}, t{3}, t{4}, t{5}, t{6}, t{7}] = mysql(cmd);
  case 8
    [t{1}, t{2}, t{3}, t{4}, t{5}, t{6}, t{7}, t{8}] = mysql(cmd);
  case 9
    [t{1}, t{2}, t{3}, t{4}, t{5}, t{6}, t{7}, t{8}, t{9}] = mysql(cmd);
  otherwise
    error('unsupported number of fields');
end

% convert the output into a structure array
for i=1:nfields
  if isnumeric(t{i})
    t{i} = num2cell(t{i});
  end
end
s = {fields{:} ; t{:}};
s = struct(s{:});
