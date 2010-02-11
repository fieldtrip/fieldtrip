function peermaster(varargin)

% PEERMASTER starts the low-level peer services and switches to
% master mode.
%
% Use as
%   peermaster(...)
%
% Optional input arguments should be specified as key-value pairs
% and can include
%   group      = string
%   hostname   = string
%
% See also PEERSLAVE, PEERRESET

% get the optional input arguments
hostname = keyval('hostname', varargin);
group    = keyval('group',    varargin);

% start the maintenance threads
warning off
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning on

if ~isempty(hostname)
  peer('hostname', hostname);
end

if ~isempty(group)
  peer('group', group);
end

% switch to master mode
peer('status', 2);

% wait some time to ensure that all current peers have been found
pause(1.5);

