function flush_data(filename, varargin)

% FLUSH_DATA removes all data from the data queue
%
% Use as
%   flush_data(filename, ...)
%
% See also FLUSH_HEADER, FLUSH_EVENT

% Copyright (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the defaults
dataformat = keyval('dataformat', varargin); if isempty(dataformat), dataformat = filetype(filename); end

switch dataformat
  case 'disp'
    % nothing to do

  case 'fcdc_buffer'
    [host, port] = filetype_check_uri(filename);
    buffer('flush_dat', [], host, port);

  case 'fcdc_mysql'
    % open the database
    [user, password, server, port] = filetype_check_uri(filename);
    if ~isempty(port)
      server = sprintf('%s:%d', server, port);
    end
    mysql('open', server, user, password);
    % remove all previous data
    cmd = 'TRUNCATE TABLE fieldtrip.data';
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
