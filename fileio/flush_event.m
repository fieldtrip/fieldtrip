function flush_event(filename, varargin)

% FLUSH_EVENT removes all events from the event queue
%
% Use as
%   flush_event(filename, ...)
%
% See also FLUSH_HEADER, FLUSH_DATA

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: flush_event.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.5  2008/10/24 08:54:13  roboos
% added support for format=matlab (i.e. simply delete the file)
%
% Revision 1.4  2008/06/19 20:50:08  roboos
% added support for fcdc_buffer
%
% Revision 1.3  2007/11/05 17:00:47  roboos
% added a comment
%
% Revision 1.2  2007/06/14 06:56:48  roboos
% created stub for flush_data and header, updated documentation
%
% Revision 1.1  2007/06/12 18:54:24  roboos
% new implementation, only mysql sofar
%

% TODO implement filtering

% set the defaults
eventformat = keyval('eventformat', varargin); if isempty(eventformat), eventformat = filetype(filename); end

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

