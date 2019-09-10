function [varargout] = istable(varargin)

% ISTABLE is a drop-in replacement for the same function that was
% introduced in MATLAB R2013b.
%
% In all MATLAB versions prior to 2013b this function returns false.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see https://github.com/fieldtrip/fieldtrip/issues/899

alternatives = which(mfilename, '-all');
if ~iscell(alternatives)
  % this is needed for octave, see https://github.com/fieldtrip/fieldtrip/pull/1171
  alternatives = {alternatives};
end

if exist(mfilename, 'builtin') || any(strncmp(alternatives, matlabroot, length(matlabroot)) & cellfun(@isempty, strfind(alternatives, fullfile('private', mfilename))))
  % remove this directory from the path
  p = fileparts(mfilename('fullpath'));
  warning('removing "%s" from your path, see http://bit.ly/2SPPjUS', p);
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
% function tf = istable(input)

tf = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {tf};
