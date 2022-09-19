function [p, d] = getcustompath

% GETCUSTOMPATH returns the directories on the MATLAB path besides those of
% the standard MATLAB toolboxes.

% Copyright (C) 2011-2012, Robert Oostenveld
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

% these are for faster processing on subsequent calls
persistent previous_path previous_argout

if isequal(path, previous_path)
  % don't do the processing again, but return the previous values from cache
  p = previous_argout{1};
  d = previous_argout{2};
  return
end

if ispc
  p = tokenize(path, ';');
else
  p = tokenize(path, ':');
end
% remove the matlab specific directories
if ispc
  s = false(size(p));
  for i=1:length(p)
    s(i) = ~strncmp(p{i}, matlabroot, length(matlabroot));
  end
else
  s = cellfun(@isempty, regexp(p, ['^' matlabroot]));
end
d = p(~s);
p = p( s);
% remove the directory containing the peer code, the worker should use its own
f = mfilename('fullpath'); % this is .../peer/private/getcustompath.m
f = fileparts(f);          % this is .../peer/private
f = fileparts(f);          % this is .../peer
p = setdiff(p, f);
% concatenate the path, using the platform specific seperator
if ispc
  p = sprintf('%s;', p{:});
  d = sprintf('%s;', d{:});
else
  p = sprintf('%s:', p{:});
  d = sprintf('%s:', d{:});
end
p = p(1:end-1); % remove the last separator
d = d(1:end-1); % remove the last separator

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_path   = path;
previous_argout = {p, d};

