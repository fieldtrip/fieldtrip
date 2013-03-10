function [ncs] = read_neuralynx_ncs(filename, begrecord, endrecord)

% READ_NEURALYNX_NCS reads a single continuous channel file
%
% Use as
%   [ncs] = read_neuralynx_ncs(filename)
%   [ncs] = read_neuralynx_ncs(filename, begrecord, endrecord)

% Copyright (C) 2005-2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
recordsize = 1044;
NRecords   = floor((ftell(fid) - headersize)/recordsize);

if NRecords>0
     
  % read out part of the dataset to detect whether there were jumps
  NRecords_to_read = min(NRecords, 100); % read out maximum 100 blocks of data
  TimeStamp        = zeros(1,NRecords_to_read,'uint64');
  ChanNumber       = zeros(1,NRecords_to_read);  
  SampFreq         = zeros(1,NRecords_to_read);
    
  for k=1:NRecords_to_read
    
    % set to the correct position    
    status = fseek(fid, headersize + (k-1)*recordsize, 'bof');        
    if status~=0
      error('cannot jump to the requested record');
    end

    % read a single continuous data record
    TimeStamp(k)    = fread(fid,   1, 'uint64=>uint64');
    ChanNumber(k)   = fread(fid,   1, 'int32');
    SampFreq(k)     = fread(fid,   1, 'int32');    
  end
  
  % for this block of data: automatically detect the gaps; 
  % there's a gap if no round off error of the sampling frequency could
  % explain the jump (which is always > one block)
  Fs       = nanmin(SampFreq);
  if Fs~=hdr.SamplingFrequency
    warning('the sampling frequency as read out from the header equals %2.2f and differs from the minimum sampling frequency as read out from the data %2.2f\n', ...
      hdr.SamplingFrequency, Fs);
  end
  
  % detect the number of timestamps per block while avoiding influencce of gaps
  d        = double(diff(TimeStamp));
  maxJump  = ceil(10^6./(Fs-1))*512;
  gapCorrectedTimeStampPerSample =  nanmean(d(d<maxJump))/512;    

  % read the timestamp from the first and last record
  if (ispc), fclose(fid); end
  ts1 = neuralynx_timestamp(filename, 1);
  tsE = neuralynx_timestamp(filename, inf);  
  if (ispc), fid = fopen(filename, 'rb', 'ieee-le'); end
  
  hdr.FirstTimeStamp  = ts1;
  hdr.LastTimeStamp   = tsE;
  
  % compare whether there's at least a block missing
  minJump = min(d);
  ts_range_predicted = (NRecords-1)*512*gapCorrectedTimeStampPerSample;
  ts_range_observed  = double(tsE-ts1);
  if abs(ts_range_predicted-ts_range_observed)>minJump
     warning('discontinuous recording, predicted number of timestamps and observed number of timestamps differ by %2.2f \n Please consult the wiki on http://fieldtrip.fcdonders.nl/getting_started/neuralynx?&#discontinuous_recordings',...
       abs(ts_range_predicted-ts_range_observed) );       
  end
      
else
  hdr.FirstTimeStamp = nan;
  hdr.LastTimeStamp  = nan;
end

if begrecord==0 && endrecord==0
  % only read the header
elseif begrecord<1
  error('cannot read before the first record');
elseif begrecord>NRecords
  error('cannot read beyond the last record')
elseif endrecord>NRecords
  endrecord = NRecords;
end

if begrecord>=1 && endrecord>=begrecord
  % rewind to the first record to be read
  status = fseek(fid, headersize + (begrecord-1)*recordsize, 'bof');
  if status~=0
    error('cannot jump to the requested record');
  end
  
  numrecord    = (endrecord-begrecord+1);
  TimeStamp    = zeros(1,numrecord,'uint64');
  ChanNumber   = zeros(1,numrecord);
  SampFreq     = zeros(1,numrecord);
  NumValidSamp = zeros(1,numrecord);
  Samp         = zeros(512,numrecord);  % this allows easy reshaping into a 1xNsamples vector
  
  for k=1:numrecord
    % read a single continuous data record
    TimeStamp(k)    = fread(fid,   1, 'uint64=>uint64');
    ChanNumber(k)   = fread(fid,   1, 'int32');
    SampFreq(k)     = fread(fid,   1, 'int32');
    NumValidSamp(k) = fread(fid,   1, 'int32');
    Samp(:,k)       = fread(fid, 512, 'int16');
    % mark the invalid samples
    Samp((NumValidSamp+1):end,k) = nan;
  end
        
  % store the record data in the output structure
  ncs.TimeStamp    = TimeStamp;
  ncs.ChanNumber   = ChanNumber;
  ncs.SampFreq     = SampFreq;
  ncs.NumValidSamp = NumValidSamp;
  % apply the scaling factor from ADBitVolts and convert to uV
  ncs.dat          = Samp * hdr.ADBitVolts * 1e6;
  
end
fclose(fid);

% store the header info in the output structure
ncs.NRecords = NRecords;
ncs.hdr      = hdr;

