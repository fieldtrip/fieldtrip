function db_open(user, password, server, port, database)

% DB_OPEN opens the connection to the database
%
% Use as
%    db_open
%    db_open(user, password, server, port, database)
%    db_open('mysql://<user>:<password>@<host>:<port>')
%
% See also DB_CLOSE, DB_SELECT, DB_INSERT, DB_SELECT_BLOB, DB_INSERT_BLOB

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

