function [dat] = read_neuralynx_ttl(filename, begsample, endsample)

% READ_NEURALYNX_TTL reads the Parallel_in values from a *.ttl file
%
% Use as
%   [dat] = read_neuralynx_ttl(filename, begsample, endsample);
%
% The *.ttl file is not a formal Neuralynx file format, but at the
% F.C. Donders Centre we use it in conjunction with Neuralynx and
% SPIKEDOWNSAMPLE.

% Copyright (C) 2006, Robert Oostenveld
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

fid = fopen(filename, 'rb', 'ieee-le');

if begsample<1
  begsample = 1;
end

if isinf(endsample)
  fseek(fid, 0, 'eof');
  endsample = (ftell(fid)-8)/4;
  fseek(fid, 0, 'bof');
end

fseek(fid,  8, 'cof');             % skip the 8-byte header with the filetype identifier
fseek(fid,  begsample-1, 'cof');   % skip to the beginning of the interesting data
dat = fread(fid, endsample-begsample+1, 'uint32=>uint32')';
if length(dat)<(endsample-begsample+1)
  error('could not read the requested data');
end

fclose(fid);

