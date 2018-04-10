function list = peerlist(status)

% PEERLIST gives information about all peers in the network, e.g. the
% number of slaves, their network configuration, etc.
%
% Use as
%   peerlist
% to get a summary of the information displayed on screen, or
%   list = peerlist
% to get all information represented in a structure.
%
% See also PEERINFO

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

% check that the  peer server threads are running
threads = true;
threads = threads && peer('announce', 'status');
threads = threads && peer('discover', 'status');
threads = threads && peer('expire',   'status');
% threads = threads && peer('tcpserver', 'status');
% threads = threads && peer('udsserver', 'status');
if ~threads
  % switch to zombie mode
  peer('status', 0);
  % start the required maintenance threads
  ws = warning('off');
  peer('announce',  'start');
  peer('discover',  'start');
  peer('expire',    'start');
  % peer('tcpserver', 'start');
  % peer('udsserver', 'start');
  warning(ws);
  % wait some time to ensure that all peers on the network have been found
  pause(1.5);
end

if nargout==0
  % this requires a complete list of all peers
  list = peer('peerlist');
  % the current field contains the job details on the busy slaves
  current = [list.current];
  if isempty(current)
    membusy = 0;
    timbusy = 0;
  else
    membusy = sum([current.memreq]);
    timbusy = sum([current.timreq]);
  end
  % give a summary
  sel = 1:numel(list);
  fprintf('there are %3d peers running in total (%d hosts, %d users)\n',length(sel), length(unique({list(sel).hostname})), length(unique({list(sel).user})));
  sel = find([list.status]==1);
  fprintf('there are %3d peers running on %2d hosts as master\n', length(sel), length(unique({list(sel).hostname})));
  sel = find([list.status]==2);
  fprintf('there are %3d peers running on %2d hosts as idle slave with %s memory available\n', length(sel), length(unique({list(sel).hostname})), print_mem(sum([list(sel).memavail])));
  sel = find([list.status]==3);
  fprintf('there are %3d peers running on %2d hosts as busy slave with %s and %s required\n', length(sel), length(unique({list(sel).hostname})), print_mem(membusy), print_tim(timbusy));
  sel = find([list.status]==0);
  fprintf('there are %3d peers running on %2d hosts as zombie\n', length(sel), length(unique({list(sel).hostname})));
end

if nargin>0
  % continue with a subset of the peers
  switch status
    case 'zombie'
      list = peer('peerlist', 0);
    case 'master'
      list = peer('peerlist', 1);
    case 'idle'
      list = peer('peerlist', 2);
    case 'busy'
      list = peer('peerlist', 3);
    otherwise
      % this requires a complete list of all peers
      list = peer('peerlist');
      % make a subset based on the host name 
      sel = zeros(size(list));
      for i=1:length(list)
        sel(i) = ~isempty(regexp(list(i).hostname, status));
      end
      list = list(find(sel));
  end
else
  % this requires a complete list of all peers
  list = peer('peerlist');
end

if nargout==0
  % the peers are listed in a random order
  % create a list which will be sorted afterward for a nice display
  strlist = cell(1,numel(list));
  for i=1:numel(list)
    switch list(i).status
      case 0
        strlist{i} = sprintf('zombie     at %s@%s:%d\n', list(i).user, list(i).hostname, list(i).port);
      case 1
        strlist{i} = sprintf('master     at %s@%s:%d\n', list(i).user, list(i).hostname, list(i).port);
      case 2
        strlist{i} = sprintf('idle slave at %s@%s:%d, memavail = %5s, timavail = %s\n', list(i).user, list(i).hostname, list(i).port, print_mem(list(i).memavail), print_tim(list(i).timavail));
      case 3
        strlist{i} = sprintf('busy slave at %s@%s:%d, working for %s, memreq = %5s, timreq = %s\n', list(i).user, list(i).hostname, list(i).port, list(i).current.user, print_mem(list(i).current.memreq), print_tim(list(i).current.timreq));
      otherwise
        error('unknown status');
    end
  end % for i
  strlist = sort(strlist);
  for i=1:numel(list)
    fprintf('%s', strlist{i});
  end
  clear list
end % if nargout

