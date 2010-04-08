function ft_flush_event(filename, varargin)

% FT_FLUSH_EVENT removes all events from the event queue
%
% Use as
%   ft_flush_event(filename, ...)
%
% See also FT_FLUSH_HEADER, FT_FLUSH_DATA

% Copyright (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% TODO implement filtering

% set the defaults
eventformat = keyval('eventformat', varargin); if isempty(eventformat), eventformat = ft_filetype(filename); end

switch eventformat
  case 'disp'
    % nothing to do

  case 'fcdc_buffer'
    [host, port] = filetype_check_uri(filename);
    buffer('flush_evt', [], host, port);

  case 'fcdc_mysql'
    % open the database
    [user, password, server, port] = filetype_check_uri(filename);
    if ~isempty(port)
      server = sprintf('%s:%d', server, port);
    end
    mysql('open', server, user, password);
    % remove all previous event information
    cmd = 'TRUNCATE TABLE fieldtrip.event';
    mysql(cmd);
    mysql('close');

  case 'matlab'
    if exist(filename, 'file')
      warning(sprintf('deleting existing file ''%s''', filename));
      delete(filename);
    end

  otherwise
    error('unsupported data format');
end

