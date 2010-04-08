function ft_flush_header(filename, varargin)

% FT_FLUSH_HEADER removes the header information from the data queue
% this also removes all data associated with the specific header.
%
% Use as
%   ft_flush_header(filename, ...)
%
% See also FT_FLUSH_DATA, FT_FLUSH_EVENT

% Copyright (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the defaults
headerformat = keyval('headerformat', varargin); if isempty(headerformat), headerformat = ft_filetype(filename); end

switch headerformat
  case 'disp'
    % nothing to do

  case 'fcdc_buffer'
    [host, port] = filetype_check_uri(filename);
    buffer('flush_hdr', [], host, port);

  case 'fcdc_mysql'
    % open the database
    [user, password, server, port] = filetype_check_uri(filename);
    if ~isempty(port)
      server = sprintf('%s:%d', server, port);
    end
    mysql('open', server, user, password);
    % remove all previous header information
    cmd = 'TRUNCATE TABLE fieldtrip.header';
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
