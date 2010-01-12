function varargout = peercellfun(fname, varargin)

% PEERCELLFUN apply a function to each element of a cell-array. The
% function execution is done in parallel on all avaialble peers.
%
% Use as
%   argout = peercellfun(fname, x1, x2, ...)
%
% This function has a number of optional arguments that have to passed
% as key-value pairs at the end of the list of input arguments. All other
% input arguments (including other key-value pairs) will be passed to the
% function to be evaluated.
%   UniformOutput  = boolean (default = false)
%
% Example
%   fname = 'power';
%   x1    = {1, 2, 3, 4, 5};
%   x2    = {2, 2, 2, 2, 2};
%   y     = peercellfun(fname, x1, x2);
%
% See also CELLFUN, PEERMASTER, PEERFEVAL, PEERGET

% locate the begin of the optional key-value arguments
optbeg = find(cellfun(@ischar, varargin));
optarg = varargin(optbeg:end);

% get the optional input arguments
UniformOutput = keyval('UniformOutput', optarg); if isempty(UniformOutput), UniformOutput = true; end
UniformOutput = istrue(UniformOutput); % convert from 'yes'/'no' into boolean value

% skip the optional key-value arguments
if ~isempty(optbeg)
  varargin = varargin(1:(optbeg-1));
end

% determine the number of input arguments and the number of jobs
numargin = numel(varargin);
numjob   = numel(varargin{1});
jobid    = nan(1, numjob);

% there are potentially errors to catch from the which() function
if isempty(which(fname))
  error('Not a valid M-file (%s).', fname);
end

% it can be difficult to determine the number of output arguments
try
  numargout = nargout(fname);
catch me
  if strcmp(me.identifier, 'MATLAB:narginout:doesNotApply')
    % e.g. in case of nargin('plus')
    numargout = 1;
  else
    rethrow(me);
  end
end

if numargout<0
  % the nargout function returns -1 in case of a variable number of output arguments
  numargout = 1;
end

% check the input arguments
for i=1:numargin
  if ~isa(varargin{i}, 'cell')
    error('input argument #%d shoudl be a cell-array', i+1);
  end
  if numel(varargin{i})~=numjob
    error('inconsistent number of elements in input #%d', i+1);
  end
end

% post all jobs
for i=1:numjob
    argin = cell(1, numargin);
    % redistribute the input arguments
    for j=1:numargin
      argin{j} = varargin{j}{i};
    end
    jobid(i) = peerfeval(fname, argin{:}, 'timeout', inf);
    clear argin
end

% get the results of all job
for indx=1:numjob
  [argout, options] = peerget(jobid(indx), 'timeout', inf, 'output', 'cell');
  % redistribute the output arguments
  for i=1:numargout
    varargout{i}{indx} = argout{i};
  end 
end

if numargout>0 && UniformOutput
  % check whether the output can be converted to a uniform one
  for i=1:numel(varargout)
    for j=1:numel(varargout{i})
      if numel(varargout{i}{j})~=1
        % this error message is consistent with the one from cellfun
        error('Non-scalar in Uniform output, at index %d, output %d. Set ''UniformOutput'' to false.', j, i);
      end
    end
  end

  % convert the output to a uniform one
  for i=1:numargout
    varargout{i} = [varargout{i}{:}];
  end
end

