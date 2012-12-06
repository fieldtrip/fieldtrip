function info = peerinfo

% PEERINFO displays information about the current peer, e.g. whether the
% maintenance threads are running, the network, group and user configuration
% and the specification of the peers that are allowed to be discovered.
%
% Use as
%   peerinfo
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
%
% $Id$
% -----------------------------------------------------------------------

info = peer('peerinfo');

if nargout==0
% display the information on screen

  % reformat the status
  switch info.status
  case 0
    status = 'zombie';
  case 1
    status = 'master';
  case 2
    status = 'idle slave';
  case 3
    status = 'busy slave';
  otherwise
    status = 'unknown';
  end

  % reformat the cell-arrays
  allowuser  = sprintf('%s, ', info.allowuser{:});
  allowgroup = sprintf('%s, ', info.allowgroup{:});
  allowhost  = sprintf('%s, ', info.allowhost{:});
  allowuser  = sprintf('{%s}', allowuser(1:end-2));
  allowgroup = sprintf('{%s}', allowgroup(1:end-2));
  allowhost  = sprintf('{%s}', allowhost(1:end-2));

  fprintf('hostid     = %u\n', info.hostid       );
  fprintf('hostname   = %s\n', info.hostname     );
  fprintf('user       = %s\n', info.user         );
  fprintf('group      = %s\n', info.group        );
  fprintf('socket     = %s\n', info.socket       );
  fprintf('port       = %d\n', info.port         );
  fprintf('status     = %s\n', status            );
  fprintf('memavail   = %u bytes\n',   info.memavail);
  fprintf('timavail   = %u seconds\n', info.timavail);
  fprintf('allowuser  = %s\n', allowuser         );
  fprintf('allowgroup = %s\n', allowgroup        );
  fprintf('allowhost  = %s\n', allowhost         );

  % give a summary of the threads
  if peer('tcpserver', 'status')
    fprintf('tcpserver thread is running\n');
  else
    fprintf('tcpserver thread is NOT running\n');
  end

  if peer('udsserver', 'status')
    fprintf('udsserver thread is running\n');
  else
    fprintf('udsserver thread is NOT running\n');
  end

  if peer('announce', 'status')
    fprintf('announce thread is running\n');
  else
    fprintf('announce thread is NOT running\n');
  end

  if peer('discover', 'status')
    fprintf('discover thread is running\n');
  else
    fprintf('discover thread is NOT running\n');
  end

  if peer('expire', 'status')
    fprintf('expire thread is running\n');
  else
    fprintf('expire thread is NOT running\n');
  end

  % give a summary of the jobs
  jobs = peer('joblist');
  list = peer('peerlist');
  fprintf('there are %d jobs in this peer''s buffer\n', length(jobs));
  for i=1:numel(jobs)
    sel = find([list.hostid] == jobs(i).hostid);
    hostname = sprintf('%s@%s:%d', list(sel).user, list(sel).hostname, list(sel).port);
    fprintf('job %d from %s with jobid %d\n', i, hostname, jobs(i).jobid);
  end

  clear info
end % if nargout

