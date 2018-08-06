function write_neuralynx_nse(filename, nse)

% WRITE_NEURALYNX_NSE writes spike timestamps and waveforms to a NSE file
% The input data should be scaled in uV.
%
% Use as
%   write_neuralynx_nse(filename, nse)
%
% See also READ_NEURALYNX_NSE

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

if ~isa(nse.TimeStamp, 'uint64')
  ft_error('timestamps should be uint64');
end

% convert the data from uV into V
nse.dat = nse.dat * 1e-6;
% scale the data and convert to 16 bits,
% this has to be done prior to writing the header
ADMaxValue = double(intmax('int16'));
ADMaxVolts = max(abs(nse.dat(:)));
if ADMaxVolts>0
  ADBitVolts = ADMaxVolts / ADMaxValue;
else
  ADBitVolts = 1;
end
nse.dat = int16(nse.dat / ADBitVolts);
% update the header with the calibration values
nse.hdr.ADBitVolts = ADBitVolts;
nse.hdr.ADMaxValue = ADMaxValue;

% construct the header
buf  = [];
buf  = [buf sprintf('######## Neuralynx Data File Header\r\n')];
buf  = [buf sprintf('## File Name: %s\r\n', filename)];
buf  = [buf sprintf('## Time Opened: (m/d/y): %s  At Time: %s\r\n', datestr(clock, 'mm/dd/yy'), datestr(clock, 'HH:MM:SS'))];
f = fieldnames(nse.hdr);
for i=1:length(f)
  v = getfield(nse.hdr, f{i});
  switch class(v)
    case 'char'
      buf = [buf sprintf('-%s\t%s\r\n', f{i}, v)];
    case 'double'
      buf = [buf sprintf('-%s\t%g\r\n', f{i}, v)];
    otherwise
      ft_error('unknown class in writing header');
  end
end

% pad the rest of the header with zeros
buf((end+1):16384) = 0;

% write the header to file
fid  = fopen(filename, 'wb', 'ieee-le');
fwrite(fid, buf);

% The format of a single electrode record is
%   int64 TimeStamp
%   int32 ScNumber
%   int32 CellNumber
%   int32 Param[0] ... Param[7]
%   int16 Samp[0] ... int16 Samp[31]

for i=1:size(nse.dat,2)
  % write a single continuous data record
  fwrite(fid, nse.TimeStamp(i) , 'int64');
  fwrite(fid, nse.ScNumber(i)  , 'int32');
  fwrite(fid, nse.CellNumber(i), 'int32');
  fwrite(fid, nse.Param(:,i)   , 'int32');
  fwrite(fid, nse.dat(:,i)     , 'int16');
end

% close the file
fclose(fid);

