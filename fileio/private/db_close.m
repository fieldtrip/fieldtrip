function db_close

% DB_CLOSE closes the connection to the database
%
% Use as
%   db_close

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: db_close.m,v $
% Revision 1.1  2007/11/07 10:50:45  roboos
% created helper functions for easy access to a MySQL database table using a structure for representing the data
%

mysql('close');

