function write_ctf_shm(msgType, msgId, sampleNumber, numSamples, numChannels, data)

% WRITE_CTF_SHM writes metainformation and data as a packet to shared memory.
% This function can be used for real-time processing of data while it is
% being acquired.
%
% Use as
%   write_ctf_shm(msgType, msgId, sampleNumber, numSamples, numChannels, data)
%
% See also READ_CTF_SHM

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
