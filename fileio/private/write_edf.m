function write_edf(filename, hdr, data)

% WRITE_EDF(filename, header, data)
%
% Writes a EDF file from the given header (only label, Fs, nChans are of interest)
% and the data (unmodified). Digital and physical limits are derived from the data
% via min and max operators. The EDF file will contain N records of 1 sample each,
% where N is the number of columns in 'data'.
%
% For sampling rates > 1 Hz, this means that the duration of one data "record"
% is less than 1s, which some EDF reading programs might complain about. At the
% same time, there is an upper limit of how big (in bytes) a record should be,
% which we could easily violate if we write the whole data as *one* record.

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

[nChans,nSamples] = size(data);
if hdr.nChans ~= nChans
  error 'Data dimension does not match header information';
end

if nChans > 9999
  error 'Cannot write more than 9999 channels to an EDF file.';
end
if nSamples > 99999999
  error 'Cannot write more than 99999999 data records (=samples) to an EDF file.';
end

labels = char(32*ones(nChans, 16));
% labels
for i=1:nChans
  ln = length(hdr.label{i});
  if ln > 16
    ft_warning('Truncating label %s to %s\n', hdr.label{i}, hdr.label{i}(1:16));
  end
  labels(i,1:ln) = hdr.label{i}(1:ln);
end

if ~isreal(data)
  ft_error('Cannot write complex-valued data');
end

scale = max(abs(data), [], 2);
scale(scale==0) = 1; % prevent division-by-zero
scale = 32767 ./ scale;

for i=1:size(data,1)
  data(i,:) = data(i,:) * scale(i);
end

data = int16(data);

% the min and max should be integers
maxV = max(data, [], 2);
minV = min(data, [], 2);

% if ~isa(data,'int16')
%   ft_warning('Warning: data type is not int16, saving to EDF might introduce round-off errors.');
%   if max(maxV) > 32767 | min(minV) < -32768
%     ft_error('Data cannot be represented as signed 16-bit integers');
%   end
%   data = int16(data);
% end

% write in data blocks of approximately one second, with an integer number of samples
blocksize = round(hdr.Fs);
nBlocks = floor(nSamples/blocksize);

digMin = sprintf('%-8i', minV);
digMax = sprintf('%-8i', maxV);

% these are tricky to print, since the '-' or 'e' messes up the width
physMin{1} = sprintf('%-8.7g', double(minV) ./ scale);
physMin{2} = sprintf('%-8.6g', double(minV) ./ scale);
physMin{3} = sprintf('%-8.5g', double(minV) ./ scale);
physMin{4} = sprintf('%-8.4g', double(minV) ./ scale);
physMin{5} = sprintf('%-8.3g', double(minV) ./ scale);
physMin{6} = sprintf('%-8.2g', double(minV) ./ scale);
physMin{7} = sprintf('%-8.1g', double(minV) ./ scale);
% take the first (most detailled) one that fits within the character space
physMin = physMin{find(cellfun(@length, physMin)==nChans*8, 1, 'first')};

physMax{1} = sprintf('%-8.7g', double(maxV) ./ scale);
physMax{2} = sprintf('%-8.6g', double(maxV) ./ scale);
physMax{3} = sprintf('%-8.5g', double(maxV) ./ scale);
physMax{4} = sprintf('%-8.4g', double(maxV) ./ scale);
physMax{5} = sprintf('%-8.3g', double(maxV) ./ scale);
physMax{6} = sprintf('%-8.2g', double(maxV) ./ scale);
physMax{7} = sprintf('%-8.1g', double(maxV) ./ scale);
% take the first (most detailled) one that fits within the character space
physMax = physMax{find(cellfun(@length, physMax)==nChans*8, 1, 'first')};

fid = fopen_or_error(filename, 'wb', 'ieee-le');
% first write fixed part
fprintf(fid, '0       ');   % version
fprintf(fid, '%-80s', '<no patient info>');
fprintf(fid, '%-80s', '<no local recording info>');

c = clock;
fprintf(fid, '%02i.%02i.%02i', c(3), c(2), mod(c(1),100)); % date as dd.mm.yy
fprintf(fid, '%02i.%02i.%02i', c(4), c(5), round(c(6))); % time as hh.mm.ss

fprintf(fid, '%-8i', 256*(1+nChans));  % number of bytes in header
fprintf(fid, '%44s', ' '); % reserved (44 spaces)
fprintf(fid, '%-8i', nBlocks); % number of data records
fprintf(fid, '%8f', blocksize/hdr.Fs); % duration of data record in seconds
fprintf(fid, '%-4i', nChans);  % number of signals = channels

fwrite(fid, labels', 'char*1'); % labels
fwrite(fid, 32*ones(80,nChans), 'uint8'); % transducer type (all spaces)
fwrite(fid, 32*ones(8,nChans), 'uint8'); % phys dimension (all spaces)
fwrite(fid, physMin', 'char*1'); % physical minimum
fwrite(fid, physMax', 'char*1'); % physical maximum
fwrite(fid, digMin', 'char*1'); % digital minimum
fwrite(fid, digMax', 'char*1'); % digital maximum
fwrite(fid, 32*ones(80,nChans), 'uint8'); % prefiltering (all spaces)
for k=1:nChans
  fprintf(fid, '%-8i', blocksize); % samples per record (each channel)
end
fwrite(fid, 32*ones(32,nChans), 'uint8'); % reserved (32 spaces / channel)

% now write data
begsample = 1;
endsample = blocksize;
while endsample<=nSamples
  fwrite(fid, data(:,begsample:endsample)', 'int16');
  begsample = begsample+blocksize;
  endsample = endsample+blocksize;
end
fclose(fid);
