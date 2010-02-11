function peermaster(varargin)

% PEERMASTER starts the low-level peer services and switches to
% master mode.
%
% See also PEERSLAVE, PEERRESET

% These need to be added
%   group      = string (default = automatic)
%   hostname   = string (default = automatic)

% start the maintenance threads
warning off
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning on

% switch to master mode
peer('status', 2);

% wait some time to ensure that all current peers have been found
pause(1.5);

