function jobid = peerfeval(varargin)

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

% the peer must be running in master mode
peer('status', 2);

% convert the input arguments into something that strmatch can work with
strargin = varargin;
strargin(~cellfun(@ischar, strargin)) = {''};

% locate the begin of the optional key-value arguments
optbeg = min([strmatch('timeout', strargin) strmatch('sleep', strargin)]);
optarg = varargin(optbeg:end);

% get the optional input arguments
timeout = keyval('timeout', optarg); if isempty(timeout), timeout = inf; end
sleep   = keyval('sleep', optarg); if isempty(sleep), sleep=0.01; end

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

% keep track of the time
stopwatch = tic;

while isempty(jobid)

  if toc(stopwatch)>timeout
    % it took too long to find a peer that was willing to execute the job
    break;
  end

  % only the peers in slave mode are interesting
  peerlist = peer('peerlist');
  peerlist = peerlist([peerlist.hoststatus]==1);

  if isempty(peerlist)
    % FIXME it would be possible in this case to execute the command locally
    error('there is no peer available as slave');
  end

  % select a random peer
  sel = floor(rand(1)*numel(peerlist))+1;
  slave = peerlist(sel);

  % there are no options to write
  options = {};

  try
    result  = peer('put', slave.hostid, varargin, options);
    jobid   = result.jobid;
  catch
    % probably the selected peer is busy, try another peer
    pause(0.01);
  end

end % while isempty(jobid)

if isempty(jobid)
  warning('none of the slave peers was willing to accept the job');
end

