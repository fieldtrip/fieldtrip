function [output, status] = ft_urlread(event_http, varargin)

% FT_URLREAD
%
% The documentation of R2016b states that urlread is not recommended.
% Use webread or webwrite instead.

method  = ft_getopt(varargin, 'method', 'get');
timeout = ft_getopt(varargin, 'timeout', 15);

% the timeout option is only available from MATLAB 2012b onward
if strcmp(method, 'get')
  if ft_platform_supports('urlread-timeout')
    [output, status] = urlread(event_http, 'timeout', timeout);
  else
    [output, status] = urlread(event_http);
  end
else
  error('method "%s" is not yet implemented', method);
end
