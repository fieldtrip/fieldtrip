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
peer('status', 2);

% convert the input arguments into something that strmatch can work with
strargin = varargin;
strargin(~cellfun(@ischar, strargin)) = {''};

% locate the begin of the optional key-value arguments
optbeg = min([strmatch('timeout', strargin) strmatch('sleep', strargin) strmatch('memreq', strargin) strmatch('cpureq', strargin) strmatch('timreq', strargin)]);
optarg = varargin(optbeg:end);

% get the optional input arguments
timeout = keyval('timeout', optarg); if isempty(timeout), timeout = inf; end
sleep   = keyval('sleep',   optarg); if isempty(sleep), sleep=0.01; end
memreq  = keyval('memreq',  optarg); if isempty(memreq), memreq=0; end
cpureq  = keyval('cpureq',  optarg); if isempty(cpureq), cpureq=0; end
timreq  = keyval('timreq',  optarg); if isempty(timreq), timreq=0; end

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
options = {'pwd', pwd, 'path', custompath};

while isempty(jobid)

  if toc(stopwatch)>timeout
    % it took too long to find a peer that was willing to execute the job
    break;
  end

  % get the full list of peers
  peerlist = peer('peerlist');

  % only peers in slave mode are interesting
  peerlist = peerlist([peerlist.hoststatus]==1);
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
  penalty    = mempenalty + 0.1* rand(1, length(peerlist));

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
z = max(x(:)) - min(x(:));
if z>0
  y = (x - min(x(:))) / z;
else
  y = (x - min(x(:)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that determines the path, excluding all Matlab toolboxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = custompath
p = tokenize(path, ':');
s = cellfun(@isempty, regexp(p, ['^' matlabroot]));
p = p(s);
p = sprintf('%s:', p{:});
p = p(1:end-1); % remove the last ':'

