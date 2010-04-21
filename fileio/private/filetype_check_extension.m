function [val] = filetype_check_extension(filename, ext)

% FILETYPE_CHECK_EXTENSION helper function to determine the file type
% by performing a case insensitive string comparison of the extension.

% Copyright (C) 2003-2006 Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if iscell(filename)
  % compare the extension of multiple files
  val = zeros(size(filename));
  for i=1:length(filename)
    val(i) = filetype_check_extension(filename{i}, ext);
  end
else
  % compare the extension of a single file
  if numel(filename)<numel(ext)
    val = 0;
  else
    val = strcmpi(filename((end-length(ext)+1):end), ext);
  end
end
return

