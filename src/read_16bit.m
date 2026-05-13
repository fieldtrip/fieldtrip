function [dat] = read_16bit(filename, offset, numwords)

% READ_16BIT read a stream of 16 bit values and converts them to doubles
% This function is designed for EDF files and is implemented as mex
% file for efficiency.
%
% Use as
%   [dat] = read_16bit(filename, offset, numwords)
%
% See also READ_24BIT

% Copyright (C) 2007, Robert Oostenveld
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

error('Could not locate the MEX file "%s.%s"', mfilename, mexext);
