function peerinfo

% PEERINFO displays information about the peers in the network and about
% the jobs that are present in this peer
%
% See also PEERLIST
%
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

list = peer('peerlist');
jobs = peer('joblist');

% display the hosts on screen, using the peerlist function
peerlist;

for i=1:numel(jobs)
  sel = find([list.hostid] == jobs(i).hostid);
  hostname = sprintf('%s@%s:%d', list(sel).user, list(sel).hostname, list(sel).hostport);
  fprintf('job from %s with jobid %d\n',  hostname, jobs(i).jobid);
end
