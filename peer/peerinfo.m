function peerinfo

% PEERINFO displays information about the peers in the network and about
% the jobs that are present in this peer
%
% See also PEERLIST

list = peer('peerlist');
jobs = peer('joblist');

% display the hosts on screen, using the peerlist function
peerlist;

for i=1:numel(jobs)
  sel = find([list.hostid] == jobs(i).hostid);
  hostname = sprintf('%s@%s:%d', list(sel).user, list(sel).hostname, list(sel).hostport);
  fprintf('job from %s with jobid %d\n',  hostname, jobs(i).jobid);
end
