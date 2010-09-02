function list = peerlist

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
  peer('status',0);     % switch to zombie mode
  pause(1.5);           % give the announce and discover some time
end

list = peer('peerlist');

if nargout==0
  % give a summary
  sel = 1:numel(list);
  fprintf('there are %3d peers running in total (%d hosts, %d users)\n',length(sel), length(unique({list(sel).hostname})), length(unique({list(sel).user}))); 
  sel = find([list.status]==1);
  fprintf('there are %3d peers running on %2d hosts as master\n',     length(sel), length(unique({list(sel).hostname})));
  sel = find([list.status]==2);
  fprintf('there are %3d peers running on %2d hosts as idle slave with %5s memory available\n', length(sel), length(unique({list(sel).hostname})), print_mem(sum([list(sel).memavail])));
  sel = find([list.status]==3);
  fprintf('there are %3d peers running on %2d hosts as busy slave with %5s memory required\n', length(sel), length(unique({list(sel).hostname})), print_mem(sum([list(sel).memavail])));
  sel = find([list.status]==0);
  fprintf('there are %3d peers running on %2d hosts as zombie\n',     length(sel), length(unique({list(sel).hostname})));

  % the peers are listed in a random order
  % create a list which will be sorted afterward for a nice display
  strlist = cell(1,numel(list));
  for i=1:numel(list)
    switch list(i).status
      case 0
        status = 'zombie     ';
      case 1
        status = 'master     ';
      case 2
        status = 'idle slave ';
      case 3
        status = 'busy slave ';
      case 4
        status = 'paused slave ';
      otherwise
        error('unknown status');
    end
    strlist{i} = sprintf('%s at %s@%s:%d, group = %s, memavail = %5s, hostid = %u\n', status, list(i).user, list(i).hostname, list(i).port, list(i).group, print_mem(list(i).memavail), list(i).hostid);
    % strlist{i} = sprintf('%s at %s@%s:%d, group = %s, memavail = %5s, timavail = %10s, hostid = %u\n', status, list(i).user, list(i).hostname, list(i).port, list(i).group, print_mem(list(i).memavail), print_tim(list(i).timavail), list(i).hostid);
  end % for i
  strlist = sort(strlist);
  for i=1:numel(list)
    fprintf('%s', strlist{i});
  end
  clear list
end

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

function str = print_tim(tim)
% partition the time in seconds into years, months, etc.
year   = 60*60*24*7*356.25;
month  = 60*60*24*356.25/12;
week   = 60*60*24*7;
day    = 60*60*24;
hour   = 60*60;
minute = 60;
org = tim; % remember the original time in seconds
y = floor(tim/year  ); tim = tim - y*year;
m = floor(tim/month ); tim = tim - m*month;
w = floor(tim/week  ); tim = tim - w*week;
d = floor(tim/day   ); tim = tim - d*day;
h = floor(tim/hour  ); tim = tim - h*hour;
m = floor(tim/minute); tim = tim - m*minute;
s = tim;
if y>=1
  str = sprintf('%.1f years', org/year);
elseif m>=1
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
  str = sprintf('%.1f seconds', org);
end

