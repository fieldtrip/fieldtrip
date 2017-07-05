function db_open(user, password, server, port, database)

% DB_OPEN opens the connection to the database
%
% Use as
%    db_open
%    db_open(user, password, server, port, database)
%    db_open('mysql://<user>:<password>@<host>:<port>')

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: db_open.m,v $
% Revision 1.2  2007/11/27 14:38:49  roboos
% added query "use database"
%
% Revision 1.1  2007/11/07 10:50:46  roboos
% created helper functions for easy access to a MySQL database table using a structure for representing the data
%

% persistent variables should be defined at the beginning of the function
persistent prev_filename

switch nargin
  case 0
    % use the default settings and combine them in a string of the form mysql://<user>:<password>@<host>:<port>
    user     = 'fieldtrip';
    password = 'fieldtrip';
    server   = 'odin';
    port     = 3306;
    filename = sprintf('mysql://<%s>:<%s>@<%s>:<%d>', user, password, server, port);
  case 1
    % the database settings are specified as mysql://<user>:<password>@<host>:<port>
    filename = user;
    [user, password, server, port] = filetype_check_uri(filename);
  case 3
    % the port was not specified, combine the settings in a string of the form mysql://<user>:<password>@<host>
    filename = sprintf('mysql://<%s>:<%s>@<%s>', user, password, server);
  case 4
    % combine the settings in a string of the form mysql://<user>:<password>@<host>:<port>
    filename = sprintf('mysql://<%s>:<%s>@<%s>:<%d>', user, password, server, port);
  case 5
    % combine the settings in a string of the form mysql://<user>:<password>@<host>:<port>
    filename = sprintf('mysql://<%s>:<%s>@<%s>:<%d>/<%s>', user, password, server, port, database);
  otherwise
    ft_error('incorrect input arguments');
end

if ~strcmp(filename, prev_filename)
  % close the database, it will be reopened furter down in the code
  mysql('close');
  prev_filename = [];
end

if ~strcmp(filename, prev_filename)
  % open the database
  if ~isempty(port)
    server = sprintf('%s:%d', server, port);
  end
  try
    mysql('open', server, user, password);
    if nargin==5
      mysql('use', database);
    end
    prev_filename = filename;
  catch
    prev_filename = [];
    ft_warning(lasterr);
  end
end

