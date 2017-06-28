function p = getcustompwd

% GETCUSTOMPWD determines the present working directory, ensuring that
% it is not the private directory of the distributed computing toolbox.

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
persistent previous_pwd previous_argout

t = pwd;

if isequal(t, previous_pwd)
  % don't do the processing again, but return the previous values from cache
  p = previous_argout;
  return
end

% don't use the present directory if it contains the peer code
% it will confuse the slave with a potentially different mex file
if strcmp(pwd, fileparts(mfilename('fullpath')))
  ft_warning('will not change directory to %s', t);
  p = [];
else
  p = t;
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_pwd    = t;
previous_argout = p;

