function [varargout] = webread(varargin)

% WEBREAD is a drop-in replacement for the function with the same
% name that was introduced in MATLAB 2014b. This function is only
% partially compatible with the original.

% Copyright (C) 2017, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

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
% function s = webread(varargin)

options = weboptions;
url = varargin{1};
s = urlread(url, 'get', options);

if isequal(s(1:2), '[{')
  % assume that it is json formatted
  split = regexp(s, '},{');
  b = [2 split+1];
  e = [split length(s)-1];
  c = cell(size(b));
  for i=1:numel(c)
    c{i} = json2struct(s(b(i):e(i)));
  end
  if numel(c)>1
    % return a cell-array with structs
    s = c;
  elseif numel(c)==1
    % return a single struct
    s = c{1};
  else
    % return the string as it is
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {s};
