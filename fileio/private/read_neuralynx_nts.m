function [nts] = read_neuralynx_nts(filename, begrecord, endrecord)

% READ_NEURALYNX_NTS reads spike timestamps
%
% Use as
%   [nts] = read_neuralynx_nts(filename)
%   [nts] = read_neuralynx_nts(filename, begrecord, endrecord)

% Copyright (C) 2006-2007, Robert Oostenveld
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

if nargin<2
  begrecord = 1;
end
if nargin<3
  endrecord = inf;
end

% the file starts with a 16*1024 bytes header in ascii, followed by a number of records
hdr = neuralynx_getheader(filename);
fid = fopen(filename, 'rb', 'ieee-le');

% determine the length of the file
fseek(fid, 0, 'eof');
headersize = 16384;
recordsize = 8;
NRecords   = floor((ftell(fid) - headersize)/recordsize);

if begrecord==0 && endrecord==0
  % only read the header  
elseif begrecord<1
  error('cannot read before the first record');
elseif begrecord>NRecords
  error('cannot read beyond the last record')
elseif endrecord>NRecords
  endrecord = NRecords;
end

if NRecords>0
  if (ispc), fclose(fid); end
  % read the timestamp from the first and last record
  hdr.FirstTimeStamp = neuralynx_timestamp(filename, 1);
  hdr.LastTimeStamp  = neuralynx_timestamp(filename, inf);
  if (ispc), fid = fopen(filename, 'rb', 'ieee-le'); end
else
  hdr.FirstTimeStamp = nan;
  hdr.LastTimeStamp  = nan;
end

if begrecord>=1 && endrecord>=begrecord
  % rewind to the first record to be read
  fseek(fid, headersize + (begrecord-1)*recordsize, 'bof');
  numrecord = (endrecord-begrecord+1);
  TimeStamp = fread(fid, numrecord, 'uint64=>uint64');

  nts.TimeStamp = TimeStamp;
end
fclose(fid);

nts.NRecords = NRecords;
nts.hdr      = hdr;
