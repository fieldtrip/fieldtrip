function peercontroller(varargin)

% PEERCONTROLLER starts the low-level peer services and switches to controller
% mode. After switching to controller mode, you can use submit jobs for
% remote execution. The server will not accept jobs to be executed,
% but does accept the output arguments of jobs that have been executed on
% other peers. Note that peercellfun will automatically execute peercontroller.
%
% Use as
%   peercontroller(...)
%
% Optional input arguments should be passed as key-value pairs. The
% following options are available to limit the peer network, i.e. to
% form sub-networks.
%   group       = string
%   allowuser   = {...}
%   allowgroup  = {...}
%   allowhost   = {...}
%   refuseuser   = {...}
%   refusegroup  = {...}
%   refusehost   = {...}
% The allow options will prevent peers that do not match the requirements
% to be added to the (dynamic) list of known peers. Consequently, these
% options limit which peers know each other. A controller will not send jobs
% to peers that it does not know. A worker will not accept jobs from a peer
% that it does not know.
%
% See also PEERWORKER, PEERRESET

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
%
% $Id$
% -----------------------------------------------------------------------

% get the optional input arguments
group       = ft_getopt(varargin, 'group');
allowuser   = ft_getopt(varargin, 'allowuser', {});
allowgroup  = ft_getopt(varargin, 'allowgroup', {});
allowhost   = ft_getopt(varargin, 'allowhost', {});
refuseuser  = ft_getopt(varargin, 'refuseuser', {});
refusegroup = ft_getopt(varargin, 'refusegroup', {});
refusehost  = ft_getopt(varargin, 'refusehost', {});

% these should be cell arrays
if ~iscell(allowuser) && ischar(allowuser)
  allowuser = {allowuser};
end
if ~iscell(allowgroup) && ischar(allowgroup)
  allowgroup = {allowgroup};
end
if ~iscell(allowhost) && ischar(allowhost)
  allowhost = {allowhost};
end
if ~iscell(refuseuser) && ischar(refuseuser)
  refuseuser = {refuseuser};
end
if ~iscell(refusegroup) && ischar(refusegroup)
  refusegroup = {refusegroup};
end
if ~iscell(refusehost) && ischar(refusehost)
  refusehost = {refusehost};
end

% this should not be user-configurable
% if ~isempty(user)
%   peer('user', user);
% end

% this should not be user-configurable
% if ~isempty(hostname)
%   peer('hostname', hostname);
% end

% the group can be specified by the user
if ~isempty(group)
  peer('group', group);
end

% switch to controller mode
peer('status', 1);

% check the current access restrictions
info   = peerinfo;
access = true;
access = access && isequal(info.allowhost,  allowhost);
access = access && isequal(info.allowuser,  allowuser);
access = access && isequal(info.allowgroup, allowgroup);
access = access && isequal(info.refusehost,  refusehost);
access = access && isequal(info.refuseuser,  refuseuser);
access = access && isequal(info.refusegroup, refusegroup);

if ~access
  % impose the updated access restrictions
  peer('allowhost',  allowhost);
  peer('allowuser',  allowuser);
  peer('allowgroup', allowgroup);
  peer('refusehost',  refusehost);
  peer('refuseuser',  refuseuser);
  peer('refusegroup', refusegroup);
end

% check the current status of the maintenance threads
threads = true;
threads = threads && peer('announce', 'status');
threads = threads && peer('discover', 'status');
threads = threads && peer('expire',   'status');
threads = threads && peer('tcpserver', 'status');
% threads = threads && peer('udsserver', 'status');

if ~threads
  % start the maintenance threads
  ws = warning('off');
  peer('announce',  'start');
  peer('discover',  'start');
  peer('expire',    'start');
  peer('tcpserver', 'start');
  % peer('udsserver', 'start');
  warning(ws);
end

if ~access || ~threads
  % wait some time to ensure that all peers on the network have been found
  pause(1.5);
end

