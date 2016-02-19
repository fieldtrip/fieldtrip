function pid = getpid()

% GETPID returns the process identifier (PID) of the current Matlab
% process.
%
% Use as
%   num = getpid;

% Copyright (C) 2011, Robert Oostenveld
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

% this is to speed up subsequent calls
persistent previous_argout
if ~isempty(previous_argout)
  pid = previous_argout;
  return
end

% this should be implemented in a MEX-file, so throw a warning here
warning('the MEX-implementation of getpid() should be used; returning a surrogate PID');

% and return a surrogate process ID
pid = round(rand(1)*1e7);

% remember for subsequent calls
previous_argout = pid;

