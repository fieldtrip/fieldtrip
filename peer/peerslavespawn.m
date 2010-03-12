function peerslavespawn(varargin)

% PEERSLAVESPAWN tries to start Matlab sessions with a peerslave
% on other computers on the network. Note that this implementation
% is site-specific. This version is tuned to the Donders Centre mentat
% cluster. You can use it as template for your own implementation.
%
% Use as
%   peerslavespawn(...)
% where the optional arguments should come in key-value pairs and 
% may include
%   maxnum       = number, don't start more than this number of slaves
%   minmemavail    = number, minimum memory requirement for slaves
%   maxmemavail    = number, maximum memory requirement for slaves
%   mintimavail    = number, minimum time requirement for slaves
%   maxtimavail    = number, maximum time requirement for slaves
%   idletim      = number, idletime the for new slaves
%
% Given the time and memory requests, this script checks
% whether a sufficient number of slaves is present and if not, it
% will try to start more slaves matching these requirements.

% get the optional input arguments
maxnum = keyval('maxnum', varargin);
minmemreq = keyval('minmemreq', varargin);
maxmemreq = keyval('maxmemreq', varargin);
mintimreq = keyval('mintimreq', varargin);
maxtimreq = keyval('maxtimreq', varargin);
idletim   = keyval('idletim', varargin);

node = textread('/opt/cluster/machines', '%s');
% node = {
%   'mentat188'
%   'mentat189'
%   'mentat190'
%   'mentat191'
%   'mentat192'
% };

% reshuffle into a random order
node = node(randperm(numel(node));

for i=1:numel(node)

  % look at the available peer-to-peer network
  list = peerlist;
  if isempty(peerlist)
    error('you can only run peerslavespawn after peermaster');
  end

  % count the number of currently suitable slaves
  suitable = true(size(list));
  suitable = status & [list.hoststatus]==1;
  suitable = status & [list.hostbusy]==0;
  suitable = status & [list.hostmemavail]>minmemavail;
  suitable = status & [list.hostmemavail]<maxmemavail;
  suitable = status & [list.hosttimavail]>mintimavail;
  suitable = status & [list.hosttimavail]<maxtimavail;
  if sum(suitable)>=maxnum
    break;
  end

  % Construct the linux command line argument which executes a shell
  % script on the remote node.  The shell script should check whether
  % the selected cluster node is suutable, and if so, start a matlab
  % session with peerslave on that node.
  cmd = sprintf('ssh public@%s peerslavespawn.sh %d %d %d %d %d', idletim, minmemavail, maxmemavail, mintimavail, maxtimavail);
  system(cmd);
end

