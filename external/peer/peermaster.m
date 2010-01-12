function peermaster

% PEERMASTER starts the low-level peer services and switches to
% master mode.
%
% See also PEERSLAVE, PEERRESET

warning off
% start the threads
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning on

% switch to master mode
peer('status', 2);

% wait some time to ensure that all current peers have been found
pause(1.5);

