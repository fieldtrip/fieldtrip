function create_buffer(port)

% CREATE_BUFFER starts the thread with the TCP server attached to the local
% Matlab instance. The TCP server will listen to the specified network
% port, and accept incoming read and write requests.
%
% Use as
%   create_buffer(port)
% where port is the TCP port to which the server listens. The default port 
% number is 1972.
% 
% See also DESTROY_BUFFER

if nargin<1
  port = 1972;
end

try
  buffer('tcpserver', 'init', 'localhost', port);
  pause(1);
catch
  if ~isempty(strfind(lasterr, 'thread is already running'))
    warning('thread is already running');
  else
    rethrow(lasterror);
  end
end
