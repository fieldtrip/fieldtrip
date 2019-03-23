function write_neuralynx_ncs(filename, ncs)

% WRITE_NEURALYNX_NCS writes continuous data to a NCS file
%
% Use as
%   write_neuralynx_ncs(filename, ncs)
%
% The input data should be scaled in uV.
%
% See also READ_NEURALYNX_NCS

% Copyright (C) 2005-2007, Robert Oostenveld
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

if ~isa(ncs.TimeStamp, 'uint64')
  ft_error('timestamps should be uint64');
end

% convert the data from uV into V
ncs.dat = ncs.dat * 1e-6;
% scale the data and convert to 16 bits,
% this has to be done prior to writing the header
ADMaxValue = double(intmax('int16'));
ADMaxVolts = max(abs(ncs.dat(:)));
if ADMaxVolts>0
  ADBitVolts = ADMaxVolts / ADMaxValue;
else
  ADBitVolts = 1;
end
ncs.dat = int16(ncs.dat / ADBitVolts);
% update the header with the calibration values
ncs.hdr.ADBitVolts = ADBitVolts;
ncs.hdr.ADMaxValue = ADMaxValue;

% construct the header
buf  = [];
buf  = [buf sprintf('######## Neuralynx Data File Header\r\n')];
buf  = [buf sprintf('## File Name: %s\r\n', filename)];
buf  = [buf sprintf('## Time Opened: (m/d/y): %s  At Time: %s\r\n', datestr(clock, 'mm/dd/yy'), datestr(clock, 'HH:MM:SS'))];
f = fieldnames(ncs.hdr);
for i=1:length(f)
  v = getfield(ncs.hdr, f{i});
  switch class(v)
    case 'char'
      buf = [buf sprintf('-%s\t%s\r\n', f{i}, v)];
    case 'double'
      buf = [buf sprintf('-%s\t%s\r\n', f{i}, num2str(v))];
    otherwise
      ft_error('unknown class in writing header');
  end
end

% pad the rest of the header with zeros
buf((end+1):16384) = 0;

% open the file and write the header
fid  = fopen_or_error(filename, 'wb', 'ieee-le');
fwrite(fid, buf);

% The format of a continuous sampled record is
%   int64 TimeStamp
%   int32 ChanNumber
%   int32 SampFreq
%   int32 NumValidSamp
%   int16 Samp[0] ... int16 Samp[511]
% Note that if NumValidSamp < 512, Samp[n], where n >= NumValidSamp, will
% contain random data.

for i=1:size(ncs.dat,2)
  % write a single continuous data record
  fwrite(fid, ncs.TimeStamp(i)   , 'uint64');
  fwrite(fid, ncs.ChanNumber(i)  , 'int32');
  fwrite(fid, ncs.SampFreq(i)    , 'int32');
  fwrite(fid, ncs.NumValidSamp(i), 'int32');
  fwrite(fid, ncs.dat(:,i)       , 'int16');
end

% close the file
fclose(fid);
