function peer

% PEER is the low-level mex file running the threads and realizing the
% communicating between the peers.
%
% peer('info')
% peer('tcpserver',   'start|stop|status')
% peer('announce',    'start|stop|status')
% peer('discover',    'start|stop|status')
% peer('expire',      'start|stop|status')
% peer('tcpport',     number)
% peer('memavail',    number)
% peer('cpuavail',    number)
% peer('timavail',    number)
% peer('clear,        jobid)
% peer('status',      0|1|2)
% peer('group',       string)
% peer('hostname',    string)
% peer('fairshare',   0|1)
% peer('allowhost',   {'host1', 'host2', ...})
% peer('allowuser',   {'user1', 'user2', ...})
% peer('allowgroup',  {'group1', 'group2', ...})
% joblist    = peer('joblist')
% peerlist   = peer('peerlist')
% peerinfo   = peer('peerinfo')
% job        = peer('put', peerid, arg, opt, ...)
% [arg, opt] = peer('get', jobid)

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

error('cannot locate mex file');

