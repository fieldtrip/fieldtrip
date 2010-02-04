function peerreset

% PEERRESET clears all jobs on the local peer server and switches to
% zombie mode.
%
% See also PEERMASTER, PEERSLAVE

peer('status', 0);

joblist = peer('joblist');
for i=1:length(joblist)
  peer('clear', joblist(i).jobid);
end

