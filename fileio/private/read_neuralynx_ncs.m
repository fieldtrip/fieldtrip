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

% Using mex files will decrease reading times 2-10-fold on big files
% I had 0.5 sec vs 3.1 sec on 20kHz 600sec file
% searching for registered precompiled Nlx2MatCSC, provided by Neuralynx
isMex = exist('Nlx2MatCSC', 'file')==3;
if isMex
    % Neuralynx mex files use C-style flags, so let's name them
    % for convinience
    READ_ALL = ones(1,5);
    flags = num2cell(diag(READ_ALL), [1,5]);
    [READ_TST READ_CHAN READ_FREQ READ_VAL READ_SAMP] = flags{:};
    % the request vector will look like: 
    %     >> TST_FLAG+FREQ_FLAG+SAMP_FLAG
    %     >> ans =
    %              1
    %              0
    %              1
    %              0
    %              1
    HEADER_NO = 0;
    HEADER_YES = 1;
    EXTRACT_RECORD_RANGE = 2;
end
% isMex = false;

% assume we have mex installed by default as it is really helpful
if ~isMex
    warning(['No Nlx2MatCSC Mex registered in the system please see '...
             'http://neuralynx.com/research_software/file_converters_and_utilities/']);
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
      
  isRead = false;
  try
      
      if ~isMex
          error('No Nlx2MatCSC Mex registered in the system');
      end
      
      % do we need ChanNumber in this section? it is not used
      [TimeStamp,  ...
       ChanNumber, ... 
       SampFreq    ] = Nlx2MatCSC(filename,...
                                  READ_TST+READ_CHAN+READ_FREQ,...
                                  HEADER_NO,... % no need in header, let FT read the header
                                  EXTRACT_RECORD_RANGE,...
                                  [1, NRecords_to_read]);
      TimeStamp = uint64(TimeStamp); % to match signature of ft_read_... output
                                     % as mex gives us doubles 
      isRead = true;
  catch Exception
      isRead = false;
%       disp(getReport(Exception));
      warning('Failed to read with Mex Nlx2MatCSC, attempting to read manually...')
      
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
      isRead = true;
  end
  
  if ~isRead 
      error('Failed to read data');
  end
  
  % explicitly sort the timestamps to deal with negative timestamp jumps that can occur
  ts1 = TimeStamp(1);
  dts = double(TimeStamp - TimeStamp(1));
  dts = unique(dts);
  dts = sort(dts);
  TimeStamp = uint64(dts) + ts1;
   
  % for this block of data: automatically detect the gaps; 
  % there's a gap if no round off error of the sampling frequency could
  % explain the jump (which is always > one block)
  Fs       = mode(double(SampFreq));
  if abs(Fs/hdr.SamplingFrequency-1)>0.01
      warning('the sampling frequency as read out from the header equals %2.2f and differs from the mode sampling frequency as read out from the data %2.2f\n', ...
      hdr.SamplingFrequency, Fs);
    
      % check which one was correct
      d = double(TimeStamp(2:end)-TimeStamp(1:end-1));  
      fsEst = 1e6./mode(d);
      indx = nearest([Fs hdr.SamplingFrequency], fsEst);
      if indx==1 
        warning('correcting the header frequency from %2.2f to %2.2f', hdr.SamplingFrequency, Fs);
        hdr.SamplingFrequency = Fs;
      end
  end
  
  % detect the number of timestamps per block while avoiding influencce of gaps
  d = double(TimeStamp(2:end)-TimeStamp(1:end-1));    
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
  % leave numrecord information here for proper synchronisation
  numrecord    = (endrecord-begrecord+1);
  
  isRead = false;
  try
      % assume we have mex installed by default as it is really helpful
      if ~isMex
          error('No Nlx2MatCSC Mex registered in the system');
      end
      
      [TimeStamp,      ...
       ChanNumber,     ... 
       SampFreq,       ...
       NumValidSamp,   ...
       Samp            ] = Nlx2MatCSC(filename,...
                                       READ_ALL,...
                                       HEADER_NO,... % no need in header
                                       EXTRACT_RECORD_RANGE,...
                                       [begrecord, endrecord]);

      TimeStamp = uint64(TimeStamp); % to match signature of ft_read_... output
      isRead = true;
  catch Exception
      isRead = false;
%       disp(getReport(Exception));
      warning('Failed to read with Mex Nlx2MatCSC, attempting to read manually...')
      
      % manual reading
      % rewind to the first record to be read
      status = fseek(fid, headersize + (begrecord-1)*recordsize, 'bof');
      if status~=0
        error('cannot jump to the requested record');
      end

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
      
      isRead = true;
  end
  
  if ~isRead 
      error('Failed to read data');
  end
  
  ts1 = TimeStamp(1);
  dts = double(TimeStamp-ts1); % no problem with doubles here as numbers are small
 
  [val,indx] = sort(dts); 
  [A,I] = unique(val); % consider only the unique values
  indx = indx(I);
    
  % store the record data in the output structure
  ncs.TimeStamp    = uint64(dts(indx)) + ts1; % convert back to original class
  ncs.ChanNumber   = ChanNumber(indx);
  ncs.SampFreq     = SampFreq(indx);
  ncs.NumValidSamp = NumValidSamp(indx);
  % apply the scaling factor from ADBitVolts and convert to uV
  ncs.dat          = Samp(:,indx) * hdr.ADBitVolts * 1e6;
end
fclose(fid);

% store the header info in the output structure
ncs.NRecords = NRecords;
ncs.hdr      = hdr;

