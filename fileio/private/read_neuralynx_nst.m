function [nst] = read_neuralynx_nst(filename, begrecord, endrecord)

% READ_NEURALYNX_NST reads a single stereotrode file
%
% Use as
%   [nst] = read_neuralynx_nst(filename)
%   [nst] = read_neuralynx_nst(filename, begrecord, endrecord)

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: read_neuralynx_nst.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.2  2008/04/29 07:52:31  roboos
% fixed windows related bug
% be consistent with begin and end timestamp in header
%
% Revision 1.1  2008/03/04 11:15:15  roboos
% new implementation, based on ntt but with only two subchannels
%

if nargin<2
  begrecord = 1;
end
if nargin<3
  endrecord = inf;
end

% The file starts with a 16*1024 bytes header in ascii, followed by a
% number of records (c.f. trials).

% The format of a tetrode record is
% int64 TimeStamp
% int32 ScNumber
% int32 CellNumber
% int32 Param[0] ƒ ƒ int32 Param[7]
% int16 ChanW[0]
% int16 ChanX[0]
% int16 ChanY[0]
% int16 ChanZ[0] ƒ ƒ
% int16 ChanW[31]
% int16 ChanX[31]
% int16 ChanY[31]
% int16 ChanZ[31]

hdr = neuralynx_getheader(filename);
fid = fopen(filename, 'rb', 'ieee-le');

% determine the length of the file
fseek(fid, 0, 'eof');
headersize = 16384;
recordsize = 304;
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
  status = fseek(fid, headersize + (begrecord-1)*recordsize, 'bof');
  if status~=0
    error('cannot jump to the requested record');
  end

  numrecord    = (endrecord-begrecord+1);
  TimeStamp    = zeros(1,numrecord,'uint64');
  ScNumber     = zeros(1,numrecord);
  CellNumber   = zeros(1,numrecord);
  Param        = zeros(8,numrecord);
  Samp         = zeros(2,32,numrecord);

  for k=1:numrecord
    TimeStamp(k)  = fread(fid, 1, 'uint64=>uint64');
    ScNumber(k)   = fread(fid, 1, 'int32');
    CellNumber(k) = fread(fid, 1, 'int32');
    Param(:,k)    = fread(fid, 8, 'int32');
    Samp(:,:,k)   = fread(fid, [2 32], 'int16'); % chan W, X
  end

  nst.TimeStamp  = TimeStamp;
  nst.ScNumber   = ScNumber;
  nst.CellNumber = CellNumber;
  nst.Param      = Param;
  % FIXME apply the scaling factor from ADBitVolts and convert to uV
  % nst.dat        = Samp * 1e6 * hdr.ADBitVolts;
  nst.dat        = Samp;
end
fclose(fid);

% store the header info in the output structure
nst.NRecords = NRecords;
nst.hdr      = hdr;
