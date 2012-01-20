function db_insert_blob(tablename, fieldname, s)

% DB_INSERT_BLOB converts a Matlab variable of arbitrary type into
% a binary stream and inserts this stream into a binary blob in the
% database table.
%
% Use as
%   db_insert_blob(tablename, fieldname, s)

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: db_insert_blob.m,v $
% Revision 1.1  2007/11/07 10:50:45  roboos
% created helper functions for easy access to a MySQL database table using a structure for representing the data
%

format = 'serialize';
cmd    = sprintf('insert into %s ( %s ) values ( "{M}" )', tablename, fieldname);
blob   = struct2msg(s, format);

% execute the query, this requires the mym version of the mysql mex file
mysql(str, blob);
