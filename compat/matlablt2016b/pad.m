function [varargout] = pad(varargin)

% PAD adds leading or trailing characters, such as spaces, to the left or
% right of an existing string.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2016b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see https://github.com/fieldtrip/fieldtrip/issues/899

if exist(mfilename, 'builtin') || any(strncmp(which(mfilename, '-all'), matlabroot, length(matlabroot)) & cellfun(@isempty, regexp(which(mfilename, '-all'), fullfile('private', mfilename))))
  % remove this directory from the path
  p = fileparts(mfilename('fullpath'));
  warning('removing %s from your path, see http://bit.ly/2SPPjUS', p);
  rmpath(p);
  % call the original MATLAB function
  if exist(mfilename, 'builtin')
    [varargout{1:nargout}] = builtin(mfilename, varargin{:});
  else
    [varargout{1:nargout}] = feval(mfilename, varargin{:});
  end
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is where the actual replacement code starts
% function s = pad(s, n, side, c)

% deal with the input arguments
if nargin==1
  [s            ] = deal(varargin{1:1});
elseif nargin==2
  [s, n         ] = deal(varargin{1:2});
elseif nargin==3
  [s, n, side   ] = deal(varargin{1:3});
elseif nargin==4
  [s, n, side, c] = deal(varargin{1:4});
else
  error('incorrect number of input arguments')
end

if nargin<4 || isempty(c)
  c = ' ';
end

if nargin<3 || isempty(side)
  side = 'right';
end

if nargin<2 || isempty(n)
  if iscell(s)
    n = max(cellfun(@length, s(:)));
  else
    n = length(s);
  end
end

if iscell(s)
  % use recursion to deal with cell-arrays
  for i=1:numel(s)
    s{i} = pad(s{i}, n, side, c);
  end

else
  % this is where the actual work happens
  assert(size(s,1)<2);
  assert(ischar(s));
  assert(numel(c)==1);
  assert(ischar(c));

  if length(s)>=n
    return
  else
    c = repmat(c, 1, n-length(s));
    switch (side)
      case 'left'
        s = [c s];
      case 'right'
        s = [s c];
      otherwise
        error('unsupported side')
    end % switch
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {s};
