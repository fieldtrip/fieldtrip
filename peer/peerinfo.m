function peerinfo

% PEERINFO displays information about the peers in the network and about
% the jobs that are present in this peer
%
% See also PEERLIST

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

% check that the required peer server threads are running
status = true;
status = status & peer('tcpserver', 'status');
status = status & peer('announce',  'status');
status = status & peer('discover',  'status');
status = status & peer('expire',    'status');
if ~status
  % start the maintenance threads
  ws = warning('off');
  peer('tcpserver', 'start');
  peer('announce',  'start');
  peer('discover',  'start');
  peer('expire',    'start');
  warning(ws)
  peer('status', 0);    % switch to zombie mode
  pause(1.5);           % give the announce and discover some time
end

list = peer('peerlist');
jobs = peer('joblist');

% give a summary
fprintf('there are %3d peers running as master\n', sum([list.hoststatus]==1));
fprintf('there are %3d peers running as idle slave\n', sum([list.hoststatus]==2));
fprintf('there are %3d peers running as busy slave\n', sum([list.hoststatus]==3));
fprintf('there are %3d peers running as zombie\n', sum([list.hoststatus]==0));

% display the hosts on screen, using the peerlist function
peerlist;

for i=1:numel(jobs)
  sel = find([list.hostid] == jobs(i).hostid);
  hostname = sprintf('%s@%s:%d', list(sel).user, list(sel).hostname, list(sel).hostport);
  fprintf('job from %s with jobid %d\n',  hostname, jobs(i).jobid);
end
