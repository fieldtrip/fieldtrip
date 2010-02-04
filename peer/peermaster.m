function peermaster(varargin)

% PEERMASTER starts the low-level peer services and switches to
% master mode.
%
% See also PEERSLAVE, PEERRESET

memavail = keyval('memavail', varargin); if isempty(memavail), memavail=intmax('uint32'); end
cpuavail = keyval('cpuavail', varargin); if isempty(cpuavail), cpuavail=intmax('uint32'); end
timavail = keyval('timavail', varargin); if isempty(timavail), timavail=intmax('uint32'); end

warning off
% start the peer server maintenance threads
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning on

% these values will be announced
peer('memavail', memavail);
peer('cpuavail', cpuavail);
peer('timavail', timavail);

% switch to master mode
peer('status', 2);

% wait some time to ensure that all current peers have been found
pause(1.5);

