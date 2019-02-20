function [varargout] = contains(varargin)

%CONTAINS True if text contains a pattern.
%   TF = contains(S,PATTERN) returns true if any element of string array S
%   contains PATTERN. TF is the same size as S.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2016b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see https://github.com/fieldtrip/fieldtrip/issues/899

if exist(mfilename, 'builtin') || strncmp(mfilename('fullpath'), matlabroot, length(matlabroot))
  % remove this directory from the path
  p = fileparts(mfilename('fullpath'));
  ft_warning('removing %s from your path, see http://bit.ly/2SPPjUS', p);
  rmpath(p);
  % call the original MATLAB function
  [varargout{1:nargout}] = builtin(mfilename, varargin{:});
  return
end

% deal with the input arguments
if nargin==1
  [s                       ] = deal(varargin{1:1});
elseif nargin==2
  [s, pattern              ] = deal(varargin{1:2});
elseif nargin==3
  [s, pattern, str         ] = deal(varargin{1:3});
elseif nargin==4
  [s, pattern, str, boolean] = deal(varargin{1:4});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is where the actual replacement code starts

if ~ischar(s) && ~iscellstr(s)
  error('the input should be either a char-array or a cell-array with chars');
end

if nargin<4
  boolean = false;
end
if nargin<3
  str = 'IgnoreCase';
end
if ~strcmpi(str, 'ignorecase')
  error('incorrect third input argument, can only be ''IgnoreCase''');
end
if ~islogical(boolean)
  error('fourth input argument should be a logical scalar');
end

if ~iscellstr(s)
  s = {s};
end

if boolean
  s       = lower(s);
  pattern = lower(pattern);
end

tf = ~cellfun(@isempty, strfind(s, pattern));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {tf};
