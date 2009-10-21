function [ncs] = read_neuralynx_ncs(filename, begrecord, endrecord)

% READ_NEURALYNX_NCS reads a single continuous channel file
%
% Use as
%   [ncs] = read_neuralynx_ncs(filename)
%   [ncs] = read_neuralynx_ncs(filename, begrecord, endrecord)

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: read_neuralynx_ncs.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.13  2008/04/29 07:52:31  roboos
% fixed windows related bug
% be consistent with begin and end timestamp in header
%
% Revision 1.12  2008/01/10 12:51:57  roboos
% ensure that it is possible to read only the header, using beg/end = 0/0
%
% Revision 1.11  2007/03/21 17:06:57  roboos
% updated the documentation
%
% Revision 1.10  2006/12/13 15:46:31  roboos
% read and keep timestamps as uint64
%
% Revision 1.9  2006/12/12 11:31:31  roboos
% cleaned up the code, made code more consistent with other neuralynx functions, moved subfunctions to seperate files, use numeric arrays instead of cell-arrays for storing the data
%
% Revision 1.8  2006/03/29 15:01:30  roboos
% fix for previous update: only apply the scaling to uV if data has been read
%
% Revision 1.7  2006/03/29 14:43:26  roboos
% scale the output data to uV, using ADBitVolts and an additional 1e6
%
% Revision 1.6  2006/03/23 18:02:13  roboos
% change endrecord from inf into the actual number present
%
% Revision 1.5  2006/03/13 16:02:55  roboos
% added first and last timestamp to the output header
%
% Revision 1.4  2006/03/10 12:32:24  roboos
% preallocate enough memory for all elements (prevent repeated memory-reallocatoin and copying)
%
% Revision 1.3  2005/11/23 08:05:23  roboos
% changed while into for-loop
% added error check for end of file in fseek
% check validity of selected begin and end record
%
% Revision 1.2  2005/09/09 12:29:38  roboos
% changed from dox to unix ascii
% changed the looping over the records
% add the number of records in the output
%
% Revision 1.1  2005/05/20 06:39:35  roboos
% new implementation of single file readers
%

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
  if (ispc), fclose(fid); end
  % read the timestamp from the first and last record
  hdr.FirstTimeStamp = neuralynx_timestamp(filename, 1);
  hdr.LastTimeStamp  = neuralynx_timestamp(filename, inf);
  if (ispc), fid = fopen(filename, 'rb', 'ieee-le'); end
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

