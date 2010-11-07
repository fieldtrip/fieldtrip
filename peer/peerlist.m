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

list = peer('peerlist');

if nargout==0
  % give a summary
  sel = 1:numel(list);
  fprintf('there are %3d peers running in total (%d hosts, %d users)\n',length(sel), length(unique({list(sel).hostname})), length(unique({list(sel).user})));
  sel = find([list.status]==1);
  fprintf('there are %3d peers running on %2d hosts as master\n', length(sel), length(unique({list(sel).hostname})));
  sel = find([list.status]==2);
  fprintf('there are %3d peers running on %2d hosts as idle slave with %s memory available\n', length(sel), length(unique({list(sel).hostname})), print_mem(sum([list(sel).memavail])));
  sel = find([list.status]==3);
  current = [list(sel).current];
  if ~isempty(current)
    memreq = sum([current.memreq]);
    timreq = sum([current.timreq]);
  else
    memreq = 0;
    timreq = 0;
  end
  fprintf('there are %3d peers running on %2d hosts as busy slave with %s and %s required\n', length(sel), length(unique({list(sel).hostname})), print_mem(memreq), print_tim(timreq));
  sel = find([list.status]==0);
  fprintf('there are %3d peers running on %2d hosts as zombie\n', length(sel), length(unique({list(sel).hostname})));
end

if nargin>0
  % continue with a subset of the peers
  switch status
    case 'zombie'
      list = list([list.status]==0);
    case 'master'
      list = list([list.status]==1);
    case 'idle'
      list = list([list.status]==2);
    case 'busy'
      list = list([list.status]==3);
    otherwise
      % this is usefull if you only want the summary
      list = [];
  end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for pretty-printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = print_mem(val)
if val<1024
  str = sprintf('%d bytes', val);
elseif val<1024^2
  str = sprintf('%.1f KB', val/1024);
elseif val<1024^3
  str = sprintf('%.1f MB', val/1024^2);
else
  str = sprintf('%.1f GB', val/1024^3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for pretty-printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = print_tim(tim)
% partition the time in seconds into years, months, etc.
year   = 60*60*24*7*356.25;
month  = 60*60*24*356.25/12;
week   = 60*60*24*7;
day    = 60*60*24;
hour   = 60*60;
minute = 60;
org = tim; % remember the original time in seconds
Y = floor(tim/year  ); tim = tim - Y*year;
M = floor(tim/month ); tim = tim - M*month;   % note capital M
w = floor(tim/week  ); tim = tim - w*week;
d = floor(tim/day   ); tim = tim - d*day;
h = floor(tim/hour  ); tim = tim - h*hour;
m = floor(tim/minute); tim = tim - m*minute;  % note small m
s = tim;
if Y>=1
  str = sprintf('%.1f years', org/year);
elseif M>=1
  str = sprintf('%.1f months', org/month);
elseif w>=1
  str = sprintf('%.1f weeks', org/week);
elseif d>=1
  str = sprintf('%.1f days', org/day);
elseif h>=1
  str = sprintf('%.1f hours', org/hour);
elseif m>=1
  str = sprintf('%.1f minutes', org/minute);
else
  % note that timreq and timavail are implemented as integers
  str = sprintf('%.0f seconds', org);
end

