function [jobid, puttime] = peerfeval(varargin)

% PEERFEVAL execute the specified function on another peer.
%
% Use as
%   jobid  = peerfeval(fname, arg1, arg2...)
%
% This function has a number of  optional arguments that have to passed
% as key-value pairs at the end of the list of input arguments. All other
% input arguments (including other key-value pairs) will be passed to the
% function to be evaluated.
%   timeout  = number, in seconds (default = inf)
%   sleep    = number, in seconds (default = 0.01)
%
% Example
%   jobid = peerfeval('unique', randn(1,10));
%   argout = peerget(jobid, 'timeout', inf);
%   disp(argout);
%
% See also FEVAL, PEERMASTER, PEERGET, PEERCELLFUN

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

% Undocumented
%   hostid   = number, only evaluate on a particular host

% check that the required peer server threads are running
status = true;
status = status & peer('tcpserver', 'status');
status = status & peer('announce',  'status');
status = status & peer('discover',  'status');
status = status & peer('expire',    'status');
if ~status
  warning('executing peermaster');
  peermaster
end

% keep track of the time
stopwatch = tic;

% the peer server must be running in master mode
peer('status', 1);

% convert the input arguments into something that strmatch can work with
strargin = varargin;
strargin(~cellfun(@ischar, strargin)) = {''};

% locate the begin of the optional key-value arguments
optbeg = false(size(strargin));
optbeg = optbeg | strcmp('timeout', strargin);
optbeg = optbeg | strcmp('sleep', strargin);
optbeg = optbeg | strcmp('memreq', strargin);
optbeg = optbeg | strcmp('cpureq', strargin);
optbeg = optbeg | strcmp('timreq', strargin);
optbeg = optbeg | strcmp('hostid', strargin);
optbeg = optbeg | strcmp('diary',  strargin);
optbeg = find(optbeg);
optarg = varargin(optbeg:end);

% get the optional input arguments
timeout = keyval('timeout', optarg); if isempty(timeout), timeout = inf; end
sleep   = keyval('sleep',   optarg); if isempty(sleep), sleep=0.01; end
memreq  = keyval('memreq',  optarg); if isempty(memreq), memreq=0; end
cpureq  = keyval('cpureq',  optarg); if isempty(cpureq), cpureq=0; end
timreq  = keyval('timreq',  optarg); if isempty(timreq), timreq=0; end
hostid  = keyval('hostid',  optarg);
diary   = keyval('diary',   optarg);

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

if isa(varargin{1}, 'function_handle')
  % convert the function handle back into a string (e.g. @plus should be 'plus')
  varargin{1} = func2str(varargin{1});
end

if isempty(which(varargin{1}))
  error('Not a valid M-file (%s).', varargin{1});
end

% start with an empty return value
jobid = [];

% pass some options that may influence remote execution
options = {'pwd', custompwd, 'path', custompath, 'diary', diary};

while isempty(jobid)
  if toc(stopwatch)>timeout
    % it took too long to find a peer that was willing to execute the job
    break;
  end
  
  % get the full list of peers
  peerlist = peer('peerlist');
  
  if ~isempty(hostid)
    % only consider peers in the user specified list
    peerlist = peerlist(ismember([peerlist.hostid], hostid));
  end
  
  % only peers in slave mode are interesting
  peerlist = peerlist([peerlist.hoststatus]==2 | [peerlist.hoststatus]==3);
  if isempty(peerlist)
    error('there is no peer available as slave');
  end
  
  % only peers with enough memory are interesting
  peerlist = peerlist([peerlist.hostmemavail] >= memreq);
  if isempty(peerlist)
    error('there are no slave peers available that meet the memory requirements');
  end
  
  % only peers with enough CPU speed are interesting
  peerlist = peerlist([peerlist.hostcpuavail] >= cpureq);
  if isempty(peerlist)
    error('there are no slave peers available that meet the CPU requirements');
  end
  
  % only peers with enough time for a single job are interesting
  peerlist = peerlist([peerlist.hosttimavail] >= timreq);
  if isempty(peerlist)
    error('there are no slave peers available to execute a job of this duration');
  end
  
  % FIXME the heuristic rule for finding the best match needs to be improved
  mempenalty = scale([peerlist.hostmemavail] - memreq);
  cpupenalty = scale([peerlist.hostcpuavail] - cpureq);
  timpenalty = scale([peerlist.hosttimavail] - timreq);
  penalty    = mempenalty + 0.1* rand(1, length(peerlist)) + ([peerlist.hoststatus]==3);
  
  % select the slave peer that has the best match with the job requirements
  [penalty, indx] = sort(penalty);
  
  % sort the peerlist according to the penalty
  peerlist = peerlist(indx);
  
  for i=1:length(peerlist)
    try
      jobid   = [];
      puttime = toc(stopwatch);
      result  = peer('put', peerlist(i).hostid, varargin, options, 'memreq', memreq, 'cpureq', cpureq, 'timreq', timreq);
      puttime = toc(stopwatch) - puttime;
      jobid   = result.jobid;
      % the peer accepted the job, there is no need to continue with the for loop
      break;
    catch
      % probably the selected peer is busy, try the next peer in line
    end
  end % for
  
  if isempty(jobid)
    % another attempt is needed, give the network some time to recover
    pause(sleep);
  end
  
end % while isempty(jobid)

if isempty(jobid)
  warning('none of the slave peers was willing to accept the job');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that scales the input values between 0 and 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = scale(x)
x = x(:);
xmin = min(x);
xmax = max(x);
if xmin==xmax
  y = (x-xmin);
else
  y = (x-xmin) / (xmax-xmin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that determines the present working directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = custompwd

% these are for faster processing on subsequent calls
persistent previous_pwd previous_argout
persistent previous_warn

if isequal(pwd, previous_pwd)
  % don't do the processing again, but return the previous values from cache
  p = previous_argout;
  return
end

% don't use the present directory if it contains the peer code
% it will confuse the slave with a potentially different mex file
if strcmp(pwd, fileparts(mfilename('fullpath')))
  if ~strcmp(previous_warn, pwd)
    warning('the peer slave will not change directory to %s', pwd);
    % warn only once
    previous_warn = pwd;
  end
  p = [];
else
  p = pwd;
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_pwd    = pwd;
previous_argout = p;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that determines the path, excluding all Matlab toolboxes
% the directories and the path on a windows computer look different than on unix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = custompath

% these are for faster processing on subsequent calls
persistent previous_path previous_argout

if isequal(path, previous_path)
  % don't do the processing again, but return the previous values from cache
  p = previous_argout;
  return
end

if ispc
  p = tokenize(path, ';');
else
  p = tokenize(path, ':');
end
% remove the matlab specific directories
if ispc
  s = false(size(p));
  for i=1:length(p)
    s(i) = ~strncmp(p{i}, matlabroot, length(matlabroot));
  end
else
  s = cellfun(@isempty, regexp(p, ['^' matlabroot]));
end
p = p(s);
% remove the directory containing the peer code, the slave should use its own
p = setdiff(p, fileparts(mfilename('fullpath')));
% concatenate the path, using the platform specific seperator
if ispc
  p = sprintf('%s;', p{:});
else
  p = sprintf('%s:', p{:});
end
p = p(1:end-1); % remove the last separator

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_path   = path;
previous_argout = p;
