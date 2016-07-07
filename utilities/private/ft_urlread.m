function [output, status] = ft_urlread(event_http)

% FT_URLREAD

% the timeout option is only available from MATLAB 2012b onward
if ft_platform_supports('urlread-timeout')
  [output, status] = urlread(event_http, 'TimeOut', 15);
else
  [output, status] = urlread(event_http);
end
