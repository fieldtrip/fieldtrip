function [varargout] = isfolder(varargin)

%ISFOLDER Determine if the input path points to a folder
%   TF = ISFOLDER(PATH) returns true if PATH points to a folder and false otherwise.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2017b.

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
% function tf = isfolder(input)

% deal with the input arguments
if nargin==1
  [input] = deal(varargin{1:1});
else
  error('incorrect number of input arguments')
end

tf = exist(input,'dir') == 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {tf};
