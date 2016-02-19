function ft_destroy_buffer

% FT_DESTROY_BUFFER stops the thread with the TCP server attached to
% the local MATLAB instance and removes all data from memory.
%
% Use as
%   ft_destroy_buffer
%
% See also FT_CREATE_BUFFER

% Copyright (C) 2010, Stefan Klanke
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

% clearing the mex file from memory will cause the function registered with
% mexAtExit to be executed. This function will then stop the threads and
% release all allocated memory to the operating system
clear buffer
