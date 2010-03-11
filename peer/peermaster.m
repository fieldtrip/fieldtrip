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
%   allowhost  = {...}
%   allowuser  = {...}
%   allowgroup = {...}
%
% See also PEERSLAVE, PEERRESET

% -----------------------------------------------------------------------
% Copyright (C) 2010, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
% -----------------------------------------------------------------------

% get the optional input arguments
hostname   = keyval('hostname',   varargin);
group      = keyval('group',      varargin);
allowhost  = keyval('allowhost',  varargin); if isempty(allowhost), allowhost = {}; end
allowuser  = keyval('allowuser',  varargin); if isempty(allowuser), allowuser = {}; end
allowgroup = keyval('allowgroup', varargin); if isempty(allowgroup), allowgroup = {}; end

% these should be cell arrays
if ~iscell(allowhost) && ischar(allowhost)
  allowhost = {allowhost};
end
if ~iscell(allowuser) && ischar(allowuser)
  allowuser = {allowuser};
end
if ~iscell(allowgroup) && ischar(allowgroup)
  allowgroup = {allowgroup};
end

% start the maintenance threads
ws = warning('off');
peer('tcpserver', 'start');
peer('announce',  'start');
peer('discover',  'start');
peer('expire',    'start');
warning(ws)

if ~isempty(hostname)
  peer('hostname', hostname);
end

if ~isempty(group)
  peer('group', group);
end

% impose access restrictions
peer('allowhost',  allowhost);
peer('allowuser',  allowuser);
peer('allowgroup', allowgroup);

% switch to master mode
peer('status', 2);

% wait some time to ensure that all current peers have been found
pause(1.5);

