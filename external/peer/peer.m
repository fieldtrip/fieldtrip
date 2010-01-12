function peer

% PEER is the low-level mex file running the threads and realizing the
% communicating between the peers.
%
% peer('info')
% peer('tcpserver', 'start|stop|status')
% peer('announce',  'start|stop|status')
% peer('discover',  'start|stop|status')
% peer('expire',    'start|stop|status')
% peer('tcpport',   number)
% peer('clear,      jobid)
% peer('status',    0|1|2)
% peer('group',     string)
% peer('allowhost',   {'host1', 'host2', ...})
% peer('allowuser',   {'user1', 'user2', ...})
% peer('allowgroup',  {'group1', 'group2', ...})
% joblist    = peer('joblist')
% peerlist   = peer('peerlist')
% job        = peer('put', peerid, arg, opt)
% job        = peer('put', peerid, arg, opt, jobid)
% [arg, opt] = peer('get', jobid)

error('cannot locate mex file');

