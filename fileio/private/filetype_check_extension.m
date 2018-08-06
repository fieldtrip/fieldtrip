function [val] = filetype_check_extension(filename, ext)

% FILETYPE_CHECK_EXTENSION helper function to determine the file type
% by performing a case insensitive string comparison of the extension.

% Copyright (C) 2003-2012 Robert Oostenveld
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

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

current_argin = {filename, ext};
if isequal(current_argin, previous_argin)
  % don't do the detection again, but return the previous value from cache
  val = previous_argout;
  return
end

if iscell(filename)
  % compare the extension of multiple files
  val = zeros(size(filename));
  for i=1:length(filename)
    val(i) = filetype_check_extension(filename{i}, ext);
  end
elseif iscell(ext)
  val = zeros(size(ext));
  for i=1:length(ext)
    val(i) = filetype_check_extension(filename, ext{i});
  end
else
  % compare the extension of a single file
  if numel(filename)<numel(ext)
    val = false;
  else
    val = strcmpi(filename((end-length(ext)+1):end), ext);
  end
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout  = val;
previous_argin  = current_argin;
previous_argout = current_argout;

return % main()
