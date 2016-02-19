function [argout] = mxDeserialize(argin)

% MXDESERIALIZE reconstructs a MATLAB object from a uint8 array suitable
% for passing down a comms channel to be reconstructed at the other end.
%
% See also MXSERIALIZE

% Copyright (C) 2005, Brad Phelan         http://xtargets.com
% Copyright (C) 2007, Robert Oostenveld   http://www.fcdonders.ru.nl
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

if verLessThan('matlab', '8.3')
  % use the original implementation of the mex file
  argout = mxDeserialize_c(argin);
else
  % use the C++ implementation of the mex file
  % see http://bugzilla.fcdonders.nl/show_bug.cgi?id=2452
  argout = mxDeserialize_cpp(argin);
end

