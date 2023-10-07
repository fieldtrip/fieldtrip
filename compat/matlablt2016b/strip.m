function [varargout] = strip(varargin)

% STRIP remove leading or trailing characters from a string.
%
% newstr = strip(str)               removes all whitespace on both sides
% newstr = strip(str, side)         removes all whitespace on the specified side
% newstr = strip(str, side, chars)  removes specified characters on the specified side
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to R2016b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to automatically remove incorrect compat folders from the path, see https://github.com/fieldtrip/fieldtrip/issues/899

if isempty(strfind(mfilename('fullpath'), matlabroot))
  % automatic cleanup does not work when FieldTrip is installed inside the MATLAB path, see https://github.com/fieldtrip/fieldtrip/issues/1527

  alternatives = which(mfilename, '-all');
  if ~iscell(alternatives)
    % this is needed for octave, see https://github.com/fieldtrip/fieldtrip/pull/1171
    alternatives = {alternatives};
  end

  keep = true(size(alternatives));
  for i=1:numel(alternatives)
    keep(i) = keep(i) && ~any(alternatives{i}=='@');  % exclude methods from classes
    keep(i) = keep(i) && alternatives{i}(end)~='p';   % exclude precompiled files
  end
  alternatives = alternatives(keep);

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

end % automatic cleanup of compat directories

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is where the actual replacement code starts
% function tf = contains(s, pattern, str, boolean)

% deal with the input arguments
if nargin==1
  str = varargin{1};
  side = 'both';
  chars = [9:13 32]; % White space characters
elseif nargin==2
  str = varargin{1};
  side = varargin{2};
  chars = [9:13 32]; % White space characters
elseif nargin==3
  str = varargin{1};
  side = varargin{2};
  chars = varargin{3};
else
  error('incorrect number of input arguments')
end

if strcmp(side, 'left') || strcmp(side, 'both')
  while any(str(1)==chars)
    str = str(2:end);
  end
end

if strcmp(side, 'right') || strcmp(side, 'both')
  while any(str(end)==chars)
    str = str(1:end-1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {str};
