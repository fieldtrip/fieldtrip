function [varargout] = split(varargin)

% SPLIT splits STR at whitespace or another user-defined character into the output
% array If STR is a string scalar or character vector, then NEWSTR is a P-by-1 string
% array or cell array where P is the number of splits.
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
% function NEWSTR = SPLIT(STR,DELIMITER)

str = varargin{1};

if nargin<2
  % split on any whitespace
  delimiter = [9:13 32];
else
  delimiter = varargin{2};
end

if ischar(str) || isstring(str)
  % find the delimiters
  sep = ismember(str, delimiter);
  % find the beginning and end of each token between the delimiters
  tokbeg = find(diff([0 ~sep])== 1);
  tokend = find(diff([~sep 0])==-1);
  newstr = cell(length(tokbeg), 1);
  for i=1:length(tokbeg)
    newstr{i} = str(tokbeg(i):tokend(i));
  end
elseif iscell()
  error('not yet implemented for cell-array input')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {newstr};
