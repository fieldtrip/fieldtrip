function setcustompwd(option_pwd)

% SETCUSTOMPWD is a helper function that updates the present working directory.

% Copyright (C) 2012, Robert Oostenveld
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

% these are for remembering the path on subsequent calls with the same input arguments
persistent previous_argin previous_pwd

if isequal(previous_argin, option_pwd) && isequal(previous_pwd, pwd)
  % no reason to change the pwd again
  return
end

if ~isempty(option_pwd)
  try
    cd(option_pwd);
  catch
    % the "catch me" syntax is broken on MATLAB74, this fixes it
    cd_error = lasterror;
    % don't throw an error, just give a warning (and hope for the best...)
    warning(cd_error.message);
  end
end

% remember the current settings for the next call
previous_argin = option_pwd;
previous_pwd   = pwd;

