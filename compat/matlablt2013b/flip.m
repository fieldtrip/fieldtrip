function [varargout] = flip(varargin)

% FLIP is a drop-in replacement for the same function that was
% introduced in MATLAB R2013b.

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
% function x = flip(x, dim)

% deal with the input arguments
if nargin==1
  [x     ] = deal(varargin{1:1});
elseif nargin==2
  [x, dim] = deal(varargin{1:2});
else
  error('incorrect number of input arguments')
end

if nargin<2 || isempty(dim)
  dim = 1;
end

n = size(x, dim);
f = n:-1:1;

switch dim
case 1
  x = x(f, :, :, :, :);
case 2
  x = x(:, f, :, :, :);
case 3
  x = x(:, :, f, :, :);
case 4
  x = x(:, :, :, f, :);
case 5
  x = x(:, :, :, :, f);
otherwise
  error('unsupported dim')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {x};
