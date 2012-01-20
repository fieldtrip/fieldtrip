function s = db_select_blob(tablename, fieldname, sel)

% DB_SELECT_BLOB selects a binary blob from a database table and converts
% it back into a Matlab variable. The variable can be of an arbitrary type.
%
% Use as
%   s = db_select_blob(tablename, fieldname)
%   s = db_select_blob(tablename, fieldname, num)
%
% The optional argument num allows you to select a specific row number.

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: db_select_blob.m,v $
% Revision 1.1  2007/11/07 10:50:46  roboos
% created helper functions for easy access to a MySQL database table using a structure for representing the data
%

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
